#include "ComponentUIFactory.h"
#include "editor/factory/DataTypeUIFactory/DataTypeFactory.h"
#include "editor/model/ModelEditorMetadata.h"
#include "editor/selection/EditorSelection.h"
#include "runtime/function/ecs/components/Model/VulkanModelComponent.h"

#include <QCheckBox>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <QFrame>
#include <QGridLayout>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QVBoxLayout>

#include <functional>
#include <iostream>

namespace {

    void ClearLayout(QLayout* layout) {
        if (layout == nullptr) {
            return;
        }

        while (QLayoutItem* item = layout->takeAt(0)) {
            if (QWidget* widget = item->widget()) {
                widget->deleteLater();
            }
            if (QLayout* childLayout = item->layout()) {
                ClearLayout(childLayout);
                delete childLayout;
            }
            delete item;
        }
    }

    void NotifyEntityDataModified() {
        Lizeral::Entity current = Lizeral::EditorSelection::Get().GetSelected();
        emit Lizeral::EditorSelection::Get().OnEntityDataModified(current);
    }

    QDoubleSpinBox* CreateFloatSpinBox(double minimum, double maximum, double step, int decimals, QWidget* parent) {
        QDoubleSpinBox* spinBox = new QDoubleSpinBox(parent);
        spinBox->setRange(minimum, maximum);
        spinBox->setSingleStep(step);
        spinBox->setDecimals(decimals);
        spinBox->setMinimumWidth(84);
        spinBox->setMinimumHeight(34);
        spinBox->setMaximumHeight(34);
        spinBox->setAlignment(Qt::AlignVCenter | Qt::AlignHCenter);
        spinBox->setButtonSymbols(QAbstractSpinBox::UpDownArrows);
        spinBox->setStyleSheet(
            "QDoubleSpinBox { padding-top: 2px; padding-bottom: 2px; }"
            "QDoubleSpinBox::up-button, QDoubleSpinBox::down-button { width: 16px; }"
        );
        return spinBox;
    }

    QDoubleSpinBox* CreateCompactFloatSpinBox(double minimum, double maximum, double step, int decimals, QWidget* parent) {
        QDoubleSpinBox* spinBox = new QDoubleSpinBox(parent);
        spinBox->setRange(minimum, maximum);
        spinBox->setSingleStep(step);
        spinBox->setDecimals(decimals);
        spinBox->setMinimumWidth(72);
        spinBox->setMinimumHeight(30);
        spinBox->setMaximumHeight(30);
        spinBox->setAlignment(Qt::AlignVCenter | Qt::AlignHCenter);
        spinBox->setButtonSymbols(QAbstractSpinBox::NoButtons);
        spinBox->setStyleSheet(
            "QDoubleSpinBox { padding: 2px 4px; }"
        );
        return spinBox;
    }

    QWidget* CreateAddressField(const QString& title, std::string& targetPath, QWidget* parent, const std::function<void()>& onEdited) {
        QWidget* widget = new QWidget(parent);
        QHBoxLayout* layout = new QHBoxLayout(widget);
        layout->setContentsMargins(0, 0, 0, 0);
        layout->setSpacing(4);

        QLabel* label = new QLabel(title, widget);
        label->setMinimumWidth(72);

        QLineEdit* pathEdit = new QLineEdit(QString::fromStdString(targetPath), widget);
        pathEdit->setPlaceholderText("Select asset path...");

        QPushButton* browseButton = new QPushButton("...", widget);
        browseButton->setFixedWidth(28);
        browseButton->setMinimumHeight(30);

        QPushButton* clearButton = new QPushButton("X", widget);
        clearButton->setFixedWidth(24);
        clearButton->setMinimumHeight(30);

        auto commitValue = [&targetPath, pathEdit, onEdited]() {
            targetPath = pathEdit->text().toStdString();
            onEdited();
        };

        QObject::connect(pathEdit, &QLineEdit::editingFinished, commitValue);
        QObject::connect(browseButton, &QPushButton::clicked, [pathEdit, commitValue, widget]() {
            const QString filePath = QFileDialog::getOpenFileName(widget, "Select Asset", pathEdit->text(), "All Files (*.*)");
            if (!filePath.isEmpty()) {
                pathEdit->setText(filePath);
                commitValue();
            }
        });
        QObject::connect(clearButton, &QPushButton::clicked, [pathEdit, commitValue]() {
            pathEdit->clear();
            commitValue();
        });

        layout->addWidget(label);
        layout->addWidget(pathEdit);
        layout->addWidget(browseButton);
        layout->addWidget(clearButton);
        return widget;
    }

    QWidget* CreateBaseColorEditor(Lizeral::VulkanMaterialSlotOverride& overrideData, QWidget* parent, const std::function<void()>& onEdited) {
        QWidget* widget = new QWidget(parent);
        QGridLayout* layout = new QGridLayout(widget);
        layout->setContentsMargins(0, 0, 0, 0);
        layout->setHorizontalSpacing(4);
        layout->setVerticalSpacing(5);

        const char* channelLabels[4] = {"R", "G", "B", "A"};
        for (int channel = 0; channel < 4; ++channel) {
            QLabel* label = new QLabel(channelLabels[channel], widget);
            label->setAlignment(Qt::AlignCenter);
            layout->addWidget(label, 0, channel);

            QDoubleSpinBox* spinBox = CreateCompactFloatSpinBox(0.0, 8.0, 0.05, 2, widget);
            spinBox->setValue(overrideData.materialInstance.factors.baseColorFactor[channel]);
            QObject::connect(spinBox, qOverload<double>(&QDoubleSpinBox::valueChanged), [channel, &overrideData, onEdited](double value) {
                overrideData.materialInstance.factors.baseColorFactor[channel] = static_cast<float>(value);
                onEdited();
            });
            layout->addWidget(spinBox, 1, channel);
        }

        return widget;
    }

    QWidget* CreateScalarOverrideRow(const QString& labelText, uint32_t maskBit, float* valuePtr, QWidget* parent, const std::function<void()>& onEdited, uint32_t* overrideMask) {
        QWidget* widget = new QWidget(parent);
        QHBoxLayout* layout = new QHBoxLayout(widget);
        layout->setContentsMargins(0, 0, 0, 0);
        layout->setSpacing(6);

        QCheckBox* enabledBox = new QCheckBox(labelText, widget);
        enabledBox->setChecked(((*overrideMask) & maskBit) != 0u);
        QObject::connect(enabledBox, &QCheckBox::checkStateChanged, [overrideMask, maskBit, onEdited](Qt::CheckState state) {
            if (state == Qt::Checked) {
                (*overrideMask) |= maskBit;
            } else {
                (*overrideMask) &= ~maskBit;
            }
            onEdited();
        });

        QDoubleSpinBox* spinBox = CreateFloatSpinBox(0.0, 1.0, 0.01, 2, widget);
        spinBox->setValue(*valuePtr);
        QObject::connect(spinBox, qOverload<double>(&QDoubleSpinBox::valueChanged), [valuePtr, onEdited](double value) {
            *valuePtr = static_cast<float>(value);
            onEdited();
        });

        layout->addWidget(enabledBox);
        layout->addStretch();
        layout->addWidget(spinBox);
        return widget;
    }

    QWidget* CreateSelectedSubMeshOverrideWidget(Lizeral::VulkanModelComponent& component, uint32_t meshAssetIndex, const Lizeral::EditorModelMetadata* metadata, QWidget* parent) {
        QGroupBox* group = new QGroupBox("Selected SubMesh", parent);
        QVBoxLayout* layout = new QVBoxLayout(group);

        const Lizeral::EditorModelMeshEntry* selectedMeshEntry = nullptr;
        if (metadata != nullptr) {
            for (const auto& meshEntry : metadata->meshEntries) {
                if (meshEntry.meshAssetIndex == meshAssetIndex) {
                    selectedMeshEntry = &meshEntry;
                    break;
                }
            }
        }

        if (selectedMeshEntry != nullptr) {
            layout->addWidget(new QLabel(QString("SubMesh: %1").arg(QString::fromStdString(selectedMeshEntry->displayName)), group));
            layout->addWidget(new QLabel(QString("Resolved Material: %1").arg(QString::fromStdString(selectedMeshEntry->materialName)), group));
        } else {
            layout->addWidget(new QLabel(QString("SubMesh Index: %1").arg(meshAssetIndex), group));
        }

        auto onEdited = [&component]() {
            component.NotifyMaterialOverridesChanged();
            NotifyEntityDataModified();
        };

        Lizeral::VulkanMaterialSlotOverride* overrideData = component.FindMaterialOverride(meshAssetIndex, 0);
        if (overrideData == nullptr) {
            QLabel* hintLabel = new QLabel("No override exists for this submesh yet.", group);
            hintLabel->setWordWrap(true);
            layout->addWidget(hintLabel);

            QPushButton* createButton = new QPushButton("Create Override", group);
            createButton->setMinimumHeight(34);
            createButton->setMinimumWidth(0);
            createButton->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
            createButton->setToolTip("Create a material override for the selected submesh.");
            QObject::connect(createButton, &QPushButton::clicked, [&component, meshAssetIndex, onEdited]() {
                component.UpsertMaterialOverride(meshAssetIndex, 0).enabled = true;
                onEdited();
                Lizeral::EditorSelection::Get().ClearSubMeshSelection();
                Lizeral::EditorSelection::Get().SelectSubMesh(Lizeral::EditorSelection::Get().GetSelected(), meshAssetIndex);
            });
            QHBoxLayout* buttonRow = new QHBoxLayout();
            buttonRow->setContentsMargins(0, 0, 0, 0);
            buttonRow->addWidget(createButton);
            layout->addLayout(buttonRow);
            return group;
        }

        QCheckBox* enabledBox = new QCheckBox("Enabled", group);
        enabledBox->setChecked(overrideData->enabled);
        QObject::connect(enabledBox, &QCheckBox::checkStateChanged, [overrideData, onEdited](Qt::CheckState state) {
            overrideData->enabled = (state == Qt::Checked);
            onEdited();
        });
        layout->addWidget(enabledBox);

        QFrame* separator = new QFrame(group);
        separator->setFrameShape(QFrame::HLine);
        layout->addWidget(separator);

        QCheckBox* baseColorMask = new QCheckBox("Override Base Color", group);
        baseColorMask->setChecked((overrideData->materialInstance.overrideMask & Lizeral::Resource::MaterialOverride_BaseColor) != 0u);
        QObject::connect(baseColorMask, &QCheckBox::checkStateChanged, [overrideData, onEdited](Qt::CheckState state) {
            if (state == Qt::Checked) {
                overrideData->materialInstance.overrideMask |= Lizeral::Resource::MaterialOverride_BaseColor;
            } else {
                overrideData->materialInstance.overrideMask &= ~Lizeral::Resource::MaterialOverride_BaseColor;
            }
            onEdited();
        });
        layout->addWidget(baseColorMask);
        layout->addWidget(CreateBaseColorEditor(*overrideData, group, onEdited));

        layout->addWidget(CreateScalarOverrideRow(
            "Override Metallic",
            Lizeral::Resource::MaterialOverride_Metallic,
            &overrideData->materialInstance.factors.metallicFactor,
            group,
            onEdited,
            &overrideData->materialInstance.overrideMask));

        layout->addWidget(CreateScalarOverrideRow(
            "Override Roughness",
            Lizeral::Resource::MaterialOverride_Roughness,
            &overrideData->materialInstance.factors.roughnessFactor,
            group,
            onEdited,
            &overrideData->materialInstance.overrideMask));

        layout->addWidget(CreateAddressField("Albedo", overrideData->textureOverrides.albedoTexturePath, group, [overrideData, onEdited]() {
            if (overrideData->textureOverrides.albedoTexturePath.empty()) {
                overrideData->materialInstance.overrideMask &= ~Lizeral::Resource::MaterialOverride_AlbedoTex;
            } else {
                overrideData->materialInstance.overrideMask |= Lizeral::Resource::MaterialOverride_AlbedoTex;
            }
            onEdited();
        }));

        layout->addWidget(CreateAddressField("Normal", overrideData->textureOverrides.normalTexturePath, group, [overrideData, onEdited]() {
            if (overrideData->textureOverrides.normalTexturePath.empty()) {
                overrideData->materialInstance.overrideMask &= ~Lizeral::Resource::MaterialOverride_NormalTex;
            } else {
                overrideData->materialInstance.overrideMask |= Lizeral::Resource::MaterialOverride_NormalTex;
            }
            onEdited();
        }));

        QLabel* ormHint = new QLabel("Packed PBR uses one ORM texture (AO / Roughness / Metallic).", group);
        ormHint->setWordWrap(true);
        layout->addWidget(ormHint);
        layout->addWidget(CreateAddressField("ORM", overrideData->textureOverrides.ormTexturePath, group, [overrideData, onEdited]() {
            if (overrideData->textureOverrides.ormTexturePath.empty()) {
                overrideData->materialInstance.overrideMask &= ~Lizeral::Resource::MaterialOverride_OrmTex;
            } else {
                overrideData->materialInstance.overrideMask |= Lizeral::Resource::MaterialOverride_OrmTex;
            }
            onEdited();
        }));

        layout->addWidget(CreateAddressField("Emissive", overrideData->textureOverrides.emissiveTexturePath, group, [overrideData, onEdited]() {
            if (overrideData->textureOverrides.emissiveTexturePath.empty()) {
                overrideData->materialInstance.overrideMask &= ~Lizeral::Resource::MaterialOverride_EmissiveTex;
            } else {
                overrideData->materialInstance.overrideMask |= Lizeral::Resource::MaterialOverride_EmissiveTex;
            }
            onEdited();
        }));

        QPushButton* removeButton = new QPushButton("Remove Override", group);
        removeButton->setMinimumHeight(34);
        QObject::connect(removeButton, &QPushButton::clicked, [&component, meshAssetIndex, onEdited]() {
            component.RemoveMaterialOverride(meshAssetIndex, 0);
            onEdited();
            Lizeral::EditorSelection::Get().ClearSubMeshSelection();
            Lizeral::EditorSelection::Get().SelectSubMesh(Lizeral::EditorSelection::Get().GetSelected(), meshAssetIndex);
        });
        layout->addWidget(removeButton);

        return group;
    }

    QWidget* CreateVulkanModelSelectionWidget(Lizeral::VulkanModelComponent& component, QWidget* parent) {
        QGroupBox* group = new QGroupBox("SubMesh Material", parent);
        QVBoxLayout* layout = new QVBoxLayout(group);

        const Lizeral::EditorModelMetadata* metadata = Lizeral::GetModelEditorMetadata(component.getModelAssetPath());
        const int32_t selectedMeshIndex = Lizeral::EditorSelection::Get().GetSelectedSubMeshIndex();
        const bool hasSelectedSubMesh = selectedMeshIndex >= 0;

        QString summary = "Select a submesh in Scene Outliner to edit its override.";
        if (metadata != nullptr && metadata->valid) {
            summary += QString("\nMeshes: %1 | Materials: %2")
                .arg(static_cast<int>(metadata->meshEntries.size()))
                .arg(static_cast<int>(metadata->materialNames.size()));
        }
        QLabel* summaryLabel = new QLabel(summary, group);
        summaryLabel->setWordWrap(true);
        layout->addWidget(summaryLabel);

        if (!hasSelectedSubMesh) {
            QLabel* hintLabel = new QLabel("No submesh is selected.", group);
            layout->addWidget(hintLabel);
            return group;
        }

        layout->addWidget(CreateSelectedSubMeshOverrideWidget(component, static_cast<uint32_t>(selectedMeshIndex), metadata, group));
        return group;
    }

} // namespace

namespace Lizeral {

    QWidget* ComponentUIFactory::CreateComponentWidget(Reflection::TypeMeta& meta, void* instance, QWidget* parent) {
        if (!instance) {
            std::cout << "[ComponentUIFactory] ERROR: Instance is null!" << std::endl;
            return nullptr;
        }

        if (!meta.isValid()) {
            std::cout << "[ComponentUIFactory] ERROR: Meta is invalid!" << std::endl;
            return nullptr;
        }

        QGroupBox* groupBox = new QGroupBox(QString::fromStdString(meta.getTypeName()), parent);
        groupBox->setStyleSheet("QGroupBox { font-weight: bold; padding-top: 15px; margin-top: 10px; }");

        QVBoxLayout* groupLayout = new QVBoxLayout(groupBox);

        auto fields = meta.getFields();
        for (auto& accessor : fields) {
            QHBoxLayout* rowLayout = new QHBoxLayout();

            QLabel* nameLabel = new QLabel(QString::fromStdString(accessor.getFieldName()), groupBox);
            nameLabel->setMinimumWidth(80);
            rowLayout->addWidget(nameLabel);

            QWidget* fieldWidget = DataTypeUIFactory::CreateFieldWidget(accessor, instance, groupBox);
            rowLayout->addWidget(fieldWidget);

            groupLayout->addLayout(rowLayout);
        }

        if (meta.getTypeName() == "VulkanModelComponent") {
            auto* modelComponent = static_cast<VulkanModelComponent*>(instance);
            groupLayout->addWidget(CreateVulkanModelSelectionWidget(*modelComponent, groupBox));
        }

        std::cout << "[ComponentUIFactory] Component widget created successfully: " << groupBox << std::endl;
        return groupBox;
    }
} // namespace Lizeral

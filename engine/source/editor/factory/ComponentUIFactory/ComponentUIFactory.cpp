#include "ComponentUIFactory.h"
#include "editor/factory/DataTypeUIFactory/DataTypeFactory.h"
#include "editor/context/EditorContext.h"
#include "editor/selection/EditorSelection.h"
#include "runtime/function/ecs/components/Model/VulkanModelComponent.h"

#include <QCheckBox>
#include <QDoubleSpinBox>
#include <QFormLayout>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QSpinBox>
#include <QVBoxLayout>

#include <functional>
#include <iostream>
#include <memory>

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
        return spinBox;
    }

    QWidget* CreateBaseColorEditor(Lizeral::VulkanMaterialSlotOverride& overrideData, QWidget* parent, const std::function<void()>& onEdited) {
        QWidget* widget = new QWidget(parent);
        QHBoxLayout* layout = new QHBoxLayout(widget);
        layout->setContentsMargins(0, 0, 0, 0);
        layout->setSpacing(4);

        for (int channel = 0; channel < 4; ++channel) {
            QDoubleSpinBox* channelSpin = CreateFloatSpinBox(0.0, 8.0, 0.05, 3, widget);
            channelSpin->setValue(overrideData.materialInstance.factors.baseColorFactor[channel]);
            QObject::connect(channelSpin, qOverload<double>(&QDoubleSpinBox::valueChanged), [channel, &overrideData, onEdited](double value) {
                overrideData.materialInstance.factors.baseColorFactor[channel] = static_cast<float>(value);
                onEdited();
            });
            layout->addWidget(channelSpin);
        }

        return widget;
    }

    QWidget* CreateTexturePathEditor(const QString& label, std::string& targetPath, QWidget* parent, const std::function<void()>& onEdited) {
        QWidget* widget = new QWidget(parent);
        QHBoxLayout* layout = new QHBoxLayout(widget);
        layout->setContentsMargins(0, 0, 0, 0);
        layout->setSpacing(4);

        QLabel* pathLabel = new QLabel(label, widget);
        pathLabel->setMinimumWidth(58);
        QLineEdit* pathEdit = new QLineEdit(QString::fromStdString(targetPath), widget);
        QObject::connect(pathEdit, &QLineEdit::editingFinished, [&targetPath, pathEdit, onEdited]() {
            targetPath = pathEdit->text().toStdString();
            onEdited();
        });

        layout->addWidget(pathLabel);
        layout->addWidget(pathEdit);
        return widget;
    }

    QWidget* CreateVulkanModelMaterialOverridesWidget(Lizeral::VulkanModelComponent& component, QWidget* parent) {
        QGroupBox* overridesGroup = new QGroupBox("Material Overrides", parent);
        QVBoxLayout* rootLayout = new QVBoxLayout(overridesGroup);

        QLabel* hintLabel = new QLabel("Per-entity overrides rebuild a dedicated material buffer the next frame.", overridesGroup);
        hintLabel->setWordWrap(true);
        rootLayout->addWidget(hintLabel);

        QWidget* overridesContainer = new QWidget(overridesGroup);
        QVBoxLayout* overridesLayout = new QVBoxLayout(overridesContainer);
        overridesLayout->setContentsMargins(0, 0, 0, 0);
        overridesLayout->setSpacing(8);
        rootLayout->addWidget(overridesContainer);

        QPushButton* addButton = new QPushButton("Add Material Override", overridesGroup);
        rootLayout->addWidget(addButton);

        auto onEdited = [&component]() {
            component.NotifyMaterialOverridesChanged();
            NotifyEntityDataModified();
        };

        std::shared_ptr<std::function<void()>> rebuild = std::make_shared<std::function<void()>>();
        *rebuild = [overridesLayout, &component, onEdited, rebuild]() {
            ClearLayout(overridesLayout);

            auto& overrides = component.GetMaterialOverrides();
            if (overrides.empty()) {
                QLabel* emptyLabel = new QLabel("No overrides yet.", nullptr);
                overridesLayout->addWidget(emptyLabel);
                return;
            }

            for (size_t index = 0; index < overrides.size(); ++index) {
                Lizeral::VulkanMaterialSlotOverride& overrideData = overrides[index];
                QGroupBox* entryGroup = new QGroupBox(QString("Override %1").arg(static_cast<int>(index)), nullptr);
                QVBoxLayout* entryLayout = new QVBoxLayout(entryGroup);

                QFormLayout* formLayout = new QFormLayout();

                QCheckBox* enabledBox = new QCheckBox(entryGroup);
                enabledBox->setChecked(overrideData.enabled);
                QObject::connect(enabledBox, &QCheckBox::checkStateChanged, [&overrideData, onEdited](Qt::CheckState state) {
                    overrideData.enabled = (state == Qt::Checked);
                    onEdited();
                });
                formLayout->addRow("Enabled", enabledBox);

                QSpinBox* meshAssetIndex = new QSpinBox(entryGroup);
                meshAssetIndex->setRange(0, 65535);
                meshAssetIndex->setValue(static_cast<int>(overrideData.meshAssetIndex));
                QObject::connect(meshAssetIndex, qOverload<int>(&QSpinBox::valueChanged), [&overrideData, onEdited](int value) {
                    overrideData.meshAssetIndex = static_cast<uint32_t>(value);
                    onEdited();
                });
                formLayout->addRow("MeshAssetIndex", meshAssetIndex);

                QSpinBox* materialSlotIndex = new QSpinBox(entryGroup);
                materialSlotIndex->setRange(0, 65535);
                materialSlotIndex->setValue(static_cast<int>(overrideData.materialSlotIndex));
                QObject::connect(materialSlotIndex, qOverload<int>(&QSpinBox::valueChanged), [&overrideData, onEdited](int value) {
                    overrideData.materialSlotIndex = static_cast<uint32_t>(value);
                    onEdited();
                });
                formLayout->addRow("MaterialSlot", materialSlotIndex);

                QLineEdit* materialAssetPath = new QLineEdit(QString::fromStdString(overrideData.materialAssetPath), entryGroup);
                QObject::connect(materialAssetPath, &QLineEdit::editingFinished, [&overrideData, materialAssetPath, onEdited]() {
                    overrideData.materialAssetPath = materialAssetPath->text().toStdString();
                    onEdited();
                });
                formLayout->addRow("MaterialAsset", materialAssetPath);

                QSpinBox* baseMaterialIndex = new QSpinBox(entryGroup);
                baseMaterialIndex->setRange(0, 65535);
                baseMaterialIndex->setValue(static_cast<int>(overrideData.materialInstance.baseMaterialIndex));
                QObject::connect(baseMaterialIndex, qOverload<int>(&QSpinBox::valueChanged), [&overrideData, onEdited](int value) {
                    overrideData.materialInstance.baseMaterialIndex = static_cast<uint32_t>(value);
                    onEdited();
                });
                formLayout->addRow("BaseMaterial", baseMaterialIndex);

                entryLayout->addLayout(formLayout);

                QWidget* baseColorWidget = CreateBaseColorEditor(overrideData, entryGroup, onEdited);
                entryLayout->addWidget(new QLabel("BaseColorFactor", entryGroup));
                entryLayout->addWidget(baseColorWidget);

                QHBoxLayout* factorsLayout = new QHBoxLayout();
                QCheckBox* baseColorMask = new QCheckBox("BaseColor", entryGroup);
                baseColorMask->setChecked((overrideData.materialInstance.overrideMask & Lizeral::Resource::MaterialOverride_BaseColor) != 0u);
                QObject::connect(baseColorMask, &QCheckBox::checkStateChanged, [&overrideData, onEdited](Qt::CheckState state) {
                    if (state == Qt::Checked) {
                        overrideData.materialInstance.overrideMask |= Lizeral::Resource::MaterialOverride_BaseColor;
                    } else {
                        overrideData.materialInstance.overrideMask &= ~Lizeral::Resource::MaterialOverride_BaseColor;
                    }
                    onEdited();
                });
                factorsLayout->addWidget(baseColorMask);

                QDoubleSpinBox* metallicSpin = CreateFloatSpinBox(0.0, 1.0, 0.01, 3, entryGroup);
                metallicSpin->setValue(overrideData.materialInstance.factors.metallicFactor);
                QObject::connect(metallicSpin, qOverload<double>(&QDoubleSpinBox::valueChanged), [&overrideData, onEdited](double value) {
                    overrideData.materialInstance.factors.metallicFactor = static_cast<float>(value);
                    onEdited();
                });
                QCheckBox* metallicMask = new QCheckBox("Metallic", entryGroup);
                metallicMask->setChecked((overrideData.materialInstance.overrideMask & Lizeral::Resource::MaterialOverride_Metallic) != 0u);
                QObject::connect(metallicMask, &QCheckBox::checkStateChanged, [&overrideData, onEdited](Qt::CheckState state) {
                    if (state == Qt::Checked) {
                        overrideData.materialInstance.overrideMask |= Lizeral::Resource::MaterialOverride_Metallic;
                    } else {
                        overrideData.materialInstance.overrideMask &= ~Lizeral::Resource::MaterialOverride_Metallic;
                    }
                    onEdited();
                });
                factorsLayout->addWidget(metallicMask);
                factorsLayout->addWidget(metallicSpin);

                QDoubleSpinBox* roughnessSpin = CreateFloatSpinBox(0.0, 1.0, 0.01, 3, entryGroup);
                roughnessSpin->setValue(overrideData.materialInstance.factors.roughnessFactor);
                QObject::connect(roughnessSpin, qOverload<double>(&QDoubleSpinBox::valueChanged), [&overrideData, onEdited](double value) {
                    overrideData.materialInstance.factors.roughnessFactor = static_cast<float>(value);
                    onEdited();
                });
                QCheckBox* roughnessMask = new QCheckBox("Roughness", entryGroup);
                roughnessMask->setChecked((overrideData.materialInstance.overrideMask & Lizeral::Resource::MaterialOverride_Roughness) != 0u);
                QObject::connect(roughnessMask, &QCheckBox::checkStateChanged, [&overrideData, onEdited](Qt::CheckState state) {
                    if (state == Qt::Checked) {
                        overrideData.materialInstance.overrideMask |= Lizeral::Resource::MaterialOverride_Roughness;
                    } else {
                        overrideData.materialInstance.overrideMask &= ~Lizeral::Resource::MaterialOverride_Roughness;
                    }
                    onEdited();
                });
                factorsLayout->addWidget(roughnessMask);
                factorsLayout->addWidget(roughnessSpin);
                entryLayout->addLayout(factorsLayout);

                auto createTextureMaskRow = [&entryLayout, &overrideData, onEdited](const QString& title, uint32_t maskBit, std::string& targetPath) {
                    QWidget* rowWidget = new QWidget(entryLayout->parentWidget());
                    QVBoxLayout* rowLayout = new QVBoxLayout(rowWidget);
                    rowLayout->setContentsMargins(0, 0, 0, 0);
                    rowLayout->setSpacing(2);

                    QCheckBox* useMask = new QCheckBox(title, rowWidget);
                    useMask->setChecked((overrideData.materialInstance.overrideMask & maskBit) != 0u);
                    QObject::connect(useMask, &QCheckBox::checkStateChanged, [&overrideData, maskBit, onEdited](Qt::CheckState state) {
                        if (state == Qt::Checked) {
                            overrideData.materialInstance.overrideMask |= maskBit;
                        } else {
                            overrideData.materialInstance.overrideMask &= ~maskBit;
                        }
                        onEdited();
                    });
                    rowLayout->addWidget(useMask);
                    rowLayout->addWidget(CreateTexturePathEditor("Path", targetPath, rowWidget, onEdited));
                    entryLayout->addWidget(rowWidget);
                };

                createTextureMaskRow("Albedo Texture", Lizeral::Resource::MaterialOverride_AlbedoTex, overrideData.textureOverrides.albedoTexturePath);
                createTextureMaskRow("Normal Texture", Lizeral::Resource::MaterialOverride_NormalTex, overrideData.textureOverrides.normalTexturePath);
                createTextureMaskRow("ORM Texture", Lizeral::Resource::MaterialOverride_OrmTex, overrideData.textureOverrides.ormTexturePath);
                createTextureMaskRow("Emissive Texture", Lizeral::Resource::MaterialOverride_EmissiveTex, overrideData.textureOverrides.emissiveTexturePath);

                QPushButton* removeButton = new QPushButton("Remove Override", entryGroup);
                QObject::connect(removeButton, &QPushButton::clicked, [&component, index, rebuild, onEdited]() {
                    auto& currentOverrides = component.GetMaterialOverrides();
                    if (index < currentOverrides.size()) {
                        currentOverrides.erase(currentOverrides.begin() + static_cast<std::ptrdiff_t>(index));
                        onEdited();
                        (*rebuild)();
                    }
                });
                entryLayout->addWidget(removeButton);

                overridesLayout->addWidget(entryGroup);
            }
        };

        QObject::connect(addButton, &QPushButton::clicked, [&component, rebuild, onEdited]() {
            uint32_t nextMeshAssetIndex = 0;
            while (component.FindMaterialOverride(nextMeshAssetIndex, 0) != nullptr) {
                ++nextMeshAssetIndex;
            }
            auto& overrideData = component.UpsertMaterialOverride(nextMeshAssetIndex, 0);
            overrideData.enabled = true;
            onEdited();
            (*rebuild)();
        });

        (*rebuild)();
        return overridesGroup;
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
            groupLayout->addWidget(CreateVulkanModelMaterialOverridesWidget(*modelComponent, groupBox));
        }

        std::cout << "[ComponentUIFactory] Component widget created successfully: " << groupBox << std::endl;
        return groupBox;
    }
} // namespace Lizeral

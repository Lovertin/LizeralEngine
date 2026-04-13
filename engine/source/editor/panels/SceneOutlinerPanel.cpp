#include "SceneOutlinerPanel.h"
#include <QVariant>
#include <QTimer>

#include "editor/selection/EditorSelection.h" 
#include "editor/model/ModelEditorMetadata.h"
#include "runtime/function/ecs/components/componentAll.h"

#include <iostream>

namespace Lizeral {

    namespace {
        constexpr int kItemTypeRole = Qt::UserRole + 1;
        constexpr int kMeshAssetIndexRole = Qt::UserRole + 2;

        enum class OutlinerItemType : int {
            Entity = 0,
            SubMesh = 1
        };
    }

    SceneOutlinerPanel::SceneOutlinerPanel(QWidget* parent) : QWidget(parent) {
        m_MainLayout = new QVBoxLayout(this);
        m_MainLayout->setContentsMargins(0, 0, 0, 0); 

        m_TreeWidget = new QTreeWidget(this);
        m_TreeWidget->setHeaderLabel("Hierarchy");
        m_TreeWidget->setContextMenuPolicy(Qt::CustomContextMenu); 
        
        m_MainLayout->addWidget(m_TreeWidget);

        // bind event
        connect(m_TreeWidget, &QTreeWidget::itemClicked, this, &SceneOutlinerPanel::OnItemClicked);
        connect(m_TreeWidget, &QTreeWidget::customContextMenuRequested, this, &SceneOutlinerPanel::ShowContextMenu);
    }

    void SceneOutlinerPanel::SetRegistry(Lizeral::Registry* registry) {
        m_Registry = registry;
        Refresh(); 
    }

    void SceneOutlinerPanel::Refresh() {
        const Entity selectedEntity = EditorSelection::Get().GetSelected();
        const int32_t selectedMeshIndex = EditorSelection::Get().GetSelectedSubMeshIndex();

        m_TreeWidget->clear();
        if (!m_Registry) return;

        auto view = m_Registry->view<NameComponent>();
        for (auto entity : view) {
            if(m_Registry->has<EditorOnlyComponent>(entity)){
                continue;
            }
            
            auto& nameComp = m_Registry->get<NameComponent>(entity);
            
            QTreeWidgetItem* item = new QTreeWidgetItem(m_TreeWidget);
            item->setText(0, QString::fromStdString(nameComp.getName()));
            item->setData(0, Qt::UserRole, QVariant::fromValue(static_cast<uint32_t>(entity)));
            item->setData(0, kItemTypeRole, static_cast<int>(OutlinerItemType::Entity));
            item->setData(0, kMeshAssetIndexRole, -1);

            if (entity == selectedEntity && selectedMeshIndex < 0) {
                m_TreeWidget->setCurrentItem(item);
            }

            if (m_Registry->has<VulkanModelComponent>(entity)) {
                auto& modelComponent = m_Registry->get<VulkanModelComponent>(entity);
                const EditorModelMetadata* metadata = GetModelEditorMetadata(modelComponent.getModelAssetPath());
                if (metadata != nullptr && metadata->valid) {
                    for (const auto& meshEntry : metadata->meshEntries) {
                        QTreeWidgetItem* childItem = new QTreeWidgetItem(item);
                        childItem->setText(0, QString::fromStdString(meshEntry.displayName));
                        childItem->setData(0, Qt::UserRole, QVariant::fromValue(static_cast<uint32_t>(entity)));
                        childItem->setData(0, kItemTypeRole, static_cast<int>(OutlinerItemType::SubMesh));
                        childItem->setData(0, kMeshAssetIndexRole, static_cast<int>(meshEntry.meshAssetIndex));

                        if (entity == selectedEntity && selectedMeshIndex == static_cast<int32_t>(meshEntry.meshAssetIndex)) {
                            item->setExpanded(true);
                            m_TreeWidget->setCurrentItem(childItem);
                        }
                    }
                }
            }
        }
    }

    void SceneOutlinerPanel::OnItemClicked(QTreeWidgetItem* item, int column) {
        if (!item) return;

        uint32_t entityId = item->data(0, Qt::UserRole).toUInt();
        Lizeral::Entity entity = static_cast<Lizeral::Entity>(entityId);
        const OutlinerItemType itemType = static_cast<OutlinerItemType>(item->data(0, kItemTypeRole).toInt());
        const int meshAssetIndex = item->data(0, kMeshAssetIndexRole).toInt();

        if (itemType == OutlinerItemType::SubMesh && meshAssetIndex >= 0) {
            EditorSelection::Get().SelectSubMesh(entity, static_cast<uint32_t>(meshAssetIndex));
        } else {
            EditorSelection::Get().SelectEntity(entity);
        }
    }

    void SceneOutlinerPanel::ShowContextMenu(const QPoint& pos) {
        if (!m_Registry) return;

        QMenu contextMenu(this);
        contextMenu.setStyleSheet("QMenu { background-color: #333333; color: white; border: 1px solid #555555; } QMenu::item:selected { background-color: #4CAF50; }");

        QAction* deleteEntityAction = nullptr;
        QTreeWidgetItem* clickedItem = m_TreeWidget->itemAt(pos);
        
        // add delete choice
        if (clickedItem) {
            deleteEntityAction = contextMenu.addAction("Delete Entity");
            contextMenu.setStyleSheet(contextMenu.styleSheet() + " QAction#DeleteAction { color: #ff6666; font-weight: bold; }");
            deleteEntityAction->setObjectName("DeleteAction");
            contextMenu.addSeparator();
        }

        QAction* createEmptyAction = contextMenu.addAction("Create Empty Entity");
        QMenu* object3DMenu = contextMenu.addMenu("3D Object");
        QAction* createCubeAction = object3DMenu->addAction("Cube");
        QMenu* lightMenu = contextMenu.addMenu("Light");
        QAction* createDirLightAction = lightMenu->addAction("Directional Light");
        QAction* createCameraAction = contextMenu.addAction("Camera");

        QAction* selectedAction = contextMenu.exec(m_TreeWidget->mapToGlobal(pos));

        if (!selectedAction) return; 

        if (selectedAction == deleteEntityAction && clickedItem) {
            uint32_t entityId = clickedItem->data(0, Qt::UserRole).toUInt();
            Lizeral::Entity targetEntity = static_cast<Lizeral::Entity>(entityId);

            QMessageBox::StandardButton reply = QMessageBox::question(this, "Delete Entity", 
                "Delete '" + clickedItem->text(0) + "' forever?", QMessageBox::Yes | QMessageBox::No);
            
            if (reply == QMessageBox::Yes) {
                // if the target Entity is delete ,refresh the Inspector 
                if (EditorSelection::Get().GetSelected() == targetEntity) {
                    EditorSelection::Get().SelectEntity(Lizeral::null_entity);
                }
                
                m_Registry->destroy(targetEntity);
                
                Refresh();
            }
        } 
        else {
            QMetaObject::invokeMethod(this, [this, selectedAction, createEmptyAction, createCubeAction, createDirLightAction, createCameraAction]() {
                
                Lizeral::Entity newEntity = m_Registry->create();
                m_Registry->emplace<TransformComponent>(newEntity); // TransformComponent is Global

                if (selectedAction == createEmptyAction) {
                    m_Registry->emplace<NameComponent>(newEntity, "Empty Entity");
                } else if (selectedAction == createCubeAction) {
                    m_Registry->emplace<NameComponent>(newEntity, "Cube");
                    m_Registry->emplace<RigidBodyComponent>(newEntity);
                    m_Registry->emplace<ColliderComponent>(newEntity);
                } else if (selectedAction == createDirLightAction) {
                    m_Registry->emplace<NameComponent>(newEntity, "Directional Light");
                    m_Registry->emplace<DirectionLightComponent>(newEntity);
                } else if (selectedAction == createCameraAction) {
                    m_Registry->emplace<NameComponent>(newEntity, "Camera");
                    m_Registry->emplace<CameraComponent>(newEntity);
                    m_Registry->emplace<CameraControlComponent>(newEntity);
                }

                Refresh();
                EditorSelection::Get().SelectEntity(newEntity);

            }, Qt::QueuedConnection);
        }
    }

} // namespace Lizeral

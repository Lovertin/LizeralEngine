#include "InspectorPanel.h"
#include "editor/factory/ComponentUIFactory/ComponentUIFactory.h"
#include "runtime/function/ecs/components/componentAll.h"
#include <iostream>

namespace Lizeral {

    InspectorPanel::InspectorPanel(QWidget* parent) : QWidget(parent) {
        // 最外层布局，永远不销毁
        QVBoxLayout* outerLayout = new QVBoxLayout(this);
        outerLayout->setContentsMargins(0, 0, 0, 0);

        // 创建第一代内容容器
        m_ContentWidget = new QWidget(this);
        m_MainLayout = new QVBoxLayout(m_ContentWidget);
        m_MainLayout->setAlignment(Qt::AlignTop); 
        
        outerLayout->addWidget(m_ContentWidget);

        InitComponentRegistry();
    }

    void InspectorPanel::InitComponentRegistry() {
        // except NameComponent as every Entity has
        RegisterComponentType<TransformComponent>("TransformComponent");
        RegisterComponentType<RigidBodyComponent>("RigidBodyComponent");
        RegisterComponentType<ColliderComponent>("ColliderComponent");
        RegisterComponentType<CameraComponent>("CameraComponent");
        RegisterComponentType<CameraControlComponent>("CameraControlComponent");
        RegisterComponentType<DirectionLightComponent>("DirectionLightComponent");
        RegisterComponentType<PointLightComponent>("PointLightComponent");
        RegisterComponentType<VulkanModelComponent>("VulkanModelComponent");
    }
    
    void InspectorPanel::ClearPanel() {
        // std::cout << "[InspectorPanel] Bulletproof Clearing Start..." << std::endl;
        if (!m_MainLayout) return;

        QLayoutItem* item;
        
        // 1. Loop through the layout, removing items from front to back
        while ((item = m_MainLayout->takeAt(0)) != nullptr) {
            // 2. If the item contains a widget
            if (QWidget* widget = item->widget()) {
                widget->setParent(nullptr); 
                widget->hide(); 
                
                widget->deleteLater(); 
            }
            // 3. We MUST delete the QLayoutItem wrapper itself.
            delete item; 
        }
        // std::cout << "[InspectorPanel] Bulletproof Clearing Finished." << std::endl;
    }

    void InspectorPanel::BindEntity(Lizeral::Entity entity) {
        if (entity == m_CurrentEntity) return;

        ClearPanel();
        m_CurrentEntity = entity;
        
        if (entity == Lizeral::null_entity || !m_Registry) return;

        auto tryDrawComponent = [&](auto* dummy, const std::string& typeName) {
            using ComponentType = std::remove_pointer_t<decltype(dummy)>;

            // find in the registry
            if (m_Registry->has<ComponentType>(entity)) {

                void* instance = &m_Registry->get<ComponentType>(entity);

                Reflection::TypeMeta meta = Reflection::TypeMeta::newMetaFromName(typeName);
                if (meta.isValid()) {
                    QWidget* compWidget = ComponentUIFactory::CreateComponentWidget(meta, instance, m_ContentWidget);
                    if (compWidget) {
                        if (typeName != "NameComponent" && typeName != "TransformComponent") {
                            QPushButton* removeBtn = new QPushButton("Remove " + QString::fromStdString(typeName));
                            removeBtn->setStyleSheet("QPushButton { background-color: #c9302c; color: white; border-radius: 3px; padding: 4px; margin: 5px; } QPushButton:hover { background-color: #d9534f; }");
                            
                            connect(removeBtn, &QPushButton::clicked, this, [this, typeName, entity]() {
                                QMessageBox::StandardButton reply = QMessageBox::question(this, "Remove Component", 
                                    "Are you sure you want to remove the " + QString::fromStdString(typeName) + "?",
                                    QMessageBox::Yes | QMessageBox::No);
                                    
                                if (reply == QMessageBox::Yes) {
                                    // erase the data
                                    m_Registry->remove<ComponentType>(entity);
                                    // refresh
                                    Lizeral::Entity target = m_CurrentEntity;
                                    m_CurrentEntity = Lizeral::null_entity;
                                    BindEntity(target);
                                }
                            });
                            
                            QVBoxLayout* groupLayout = static_cast<QVBoxLayout*>(compWidget->layout());
                            if (groupLayout) groupLayout->addWidget(removeBtn);
                        }
                        m_MainLayout->addWidget(compWidget);
                    }
                }
            }
        };

        // ordered
        tryDrawComponent((NameComponent*)nullptr, "NameComponent");
        tryDrawComponent((TransformComponent*)nullptr, "TransformComponent");
        tryDrawComponent((RigidBodyComponent*)nullptr, "RigidBodyComponent");
        tryDrawComponent((ColliderComponent*)nullptr,  "ColliderComponent");
        tryDrawComponent((CameraComponent*)nullptr,    "CameraComponent");
        tryDrawComponent((CameraControlComponent*)nullptr, "CameraControlComponent");
        tryDrawComponent((DirectionLightComponent*)nullptr, "DirectionLightComponent");
        tryDrawComponent((PointLightComponent*)nullptr, "PointLightComponent");
        tryDrawComponent((VulkanModelComponent*)nullptr, "VulkanModelComponent");

        QPushButton* addCompBtn = new QPushButton("Add Component", m_ContentWidget);
        addCompBtn->setStyleSheet("QPushButton { margin-top: 15px; padding: 8px; font-weight: bold; background-color: #4CAF50; color: white; border-radius: 4px; } QPushButton:hover { background-color: #45a049; }");
        
        connect(addCompBtn, &QPushButton::clicked, this, [this, addCompBtn]() {
            ShowAddComponentMenu(addCompBtn);
        });

        m_MainLayout->addWidget(addCompBtn);

        m_MainLayout->addStretch(); 
    }


    void InspectorPanel::ShowAddComponentMenu(QPushButton* button) {
        QMenu menu(nullptr);
        menu.setStyleSheet("QMenu { background-color: #333333; color: white; border: 1px solid #555555; } QMenu::item:selected { background-color: #4CAF50; }");

        bool hasAnyAvailable = false;

        for (const auto& compAction : m_AvailableComponents) {
            if (!compAction.hasComponent(m_Registry, m_CurrentEntity)) {
                hasAnyAvailable = true;
                QAction* action = menu.addAction(QString::fromStdString(compAction.name));

                connect(action, &QAction::triggered, this, [this, compAction]() {
                    compAction.addComponent(m_Registry, m_CurrentEntity);

                    Lizeral::Entity target = m_CurrentEntity;
                    m_CurrentEntity = Lizeral::null_entity; 
                    BindEntity(target);
                });
            }
        }

        if (!hasAnyAvailable) {
            QAction* emptyAction = menu.addAction("No components to add");
            emptyAction->setEnabled(false);
        }

        menu.exec(button->mapToGlobal(QPoint(0, button->height())));
    }

} // namespace Lizeral

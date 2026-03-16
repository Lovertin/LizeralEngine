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
        // 这里不要注册 NameComponent，因为每个实体默认都必须有名字
        RegisterComponentType<TransformComponent>("TransformComponent");
        RegisterComponentType<RigidBodyComponent>("RigidBodyComponent");
        RegisterComponentType<ColliderComponent>("ColliderComponent");
        RegisterComponentType<CameraComponent>("CameraComponent");
        RegisterComponentType<CameraControlComponent>("CameraControlComponent");
        RegisterComponentType<DirectionLightComponent>("DirectionLightComponent");
        RegisterComponentType<ModelComponent>("ModelComponent");
    }
    
    void InspectorPanel::ClearPanel() {
        // std::cout << "[InspectorPanel] Bulletproof Clearing Start..." << std::endl;
        if (!m_MainLayout) return;

        QLayoutItem* item;
        
        // 1. Loop through the layout, removing items from front to back
        while ((item = m_MainLayout->takeAt(0)) != nullptr) {
            // 2. If the item contains a widget
            if (QWidget* widget = item->widget()) {
                // 彻底断开与界面的连接
                widget->setParent(nullptr); 
                widget->hide(); 
                
                // 延迟到事件循环清理，极其安全
                widget->deleteLater(); 
            }
            // 3. We MUST delete the QLayoutItem wrapper itself.
            delete item; 
        }
        // std::cout << "[InspectorPanel] Bulletproof Clearing Finished." << std::endl;
    }

    void InspectorPanel::BindEntity(Lizeral::Entity entity) {
        // 防止重复绑定导致无意义的重绘
        if (entity == m_CurrentEntity) return;

        ClearPanel();
        m_CurrentEntity = entity;
        
        if (entity == Lizeral::null_entity || !m_Registry) return;

        // ==========================================================
        // 核心魔法：查询 ECS 并动态桥接反射系统
        // ==========================================================
        auto tryDrawComponent = [&](auto* dummy, const std::string& typeName) {
            using ComponentType = std::remove_pointer_t<decltype(dummy)>;

            // 1. 去 ECS 里查，这个实体有没有这个组件？
            if (m_Registry->has<ComponentType>(entity)) {
                // 2. 拿到真实的内存地址
                void* instance = &m_Registry->get<ComponentType>(entity);
                
                // 3. 拿到反射元数据
                Reflection::TypeMeta meta = Reflection::TypeMeta::newMetaFromName(typeName);
                if (meta.isValid()) {
                    // 4. 丢给组件工厂去画大框框，注意父节点要传 m_ContentWidget
                    QWidget* compWidget = ComponentUIFactory::CreateComponentWidget(meta, instance, m_ContentWidget);
                    if (compWidget) {
                        if (typeName != "NameComponent" && typeName != "TransformComponent") {
                            QPushButton* removeBtn = new QPushButton("Remove " + QString::fromStdString(typeName));
                            removeBtn->setStyleSheet("QPushButton { background-color: #c9302c; color: white; border-radius: 3px; padding: 4px; margin: 5px; } QPushButton:hover { background-color: #d9534f; }");
                            
                            connect(removeBtn, &QPushButton::clicked, this, [this, typeName, entity]() {
                                // 弹窗二次确认
                                QMessageBox::StandardButton reply = QMessageBox::question(this, "Remove Component", 
                                    "Are you sure you want to remove the " + QString::fromStdString(typeName) + "?",
                                    QMessageBox::Yes | QMessageBox::No);
                                    
                                if (reply == QMessageBox::Yes) {
                                    // 1. 底层擦除数据
                                    m_Registry->remove<ComponentType>(entity);
                                    // 2. 刷新面板
                                    Lizeral::Entity target = m_CurrentEntity;
                                    m_CurrentEntity = Lizeral::null_entity;
                                    BindEntity(target);
                                }
                            });
                            
                            // 强行把按钮塞进 ComponentUIFactory 生成的 GroupBox 布局的最下面
                            QVBoxLayout* groupLayout = static_cast<QVBoxLayout*>(compWidget->layout());
                            if (groupLayout) groupLayout->addWidget(removeBtn);
                        }
                        m_MainLayout->addWidget(compWidget);
                    }
                }
            }
        };

        // 依次排查并绘制所有已知的组件
        // （这里的顺序，就是 Inspector 面板里自上而下的显示顺序）
        tryDrawComponent((NameComponent*)nullptr, "NameComponent");
        tryDrawComponent((TransformComponent*)nullptr, "TransformComponent");
        tryDrawComponent((RigidBodyComponent*)nullptr, "RigidBodyComponent");
        tryDrawComponent((ColliderComponent*)nullptr,  "ColliderComponent");
        tryDrawComponent((CameraComponent*)nullptr,    "CameraComponent");
        tryDrawComponent((CameraControlComponent*)nullptr, "CameraControlComponent");
        tryDrawComponent((DirectionLightComponent*)nullptr, "DirectionLightComponent");
        tryDrawComponent((VulkanModelComponent*)nullptr, "ModelComponent");
        // tryDrawComponent((ModelComponent*)nullptr,     "ModelComponent");

        QPushButton* addCompBtn = new QPushButton("Add Component", m_ContentWidget);
        addCompBtn->setStyleSheet("QPushButton { margin-top: 15px; padding: 8px; font-weight: bold; background-color: #4CAF50; color: white; border-radius: 4px; } QPushButton:hover { background-color: #45a049; }");
        
        connect(addCompBtn, &QPushButton::clicked, this, [this, addCompBtn]() {
            ShowAddComponentMenu(addCompBtn);
        });

        m_MainLayout->addWidget(addCompBtn);

        // 在最底部加个弹簧，把所有组件往上顶
        m_MainLayout->addStretch(); 
    }


    void InspectorPanel::ShowAddComponentMenu(QPushButton* button) {
        QMenu menu(this);
        menu.setStyleSheet("QMenu { background-color: #333333; color: white; border: 1px solid #555555; } QMenu::item:selected { background-color: #4CAF50; }");

        bool hasAnyAvailable = false;

        // 遍历我们注册的所有组件
        for (const auto& compAction : m_AvailableComponents) {
            // 【核心逻辑】：如果这个实体还没有这个组件，才把它加到菜单里！
            if (!compAction.hasComponent(m_Registry, m_CurrentEntity)) {
                hasAnyAvailable = true;
                QAction* action = menu.addAction(QString::fromStdString(compAction.name));
                
                // 当用户点击菜单项时触发
                connect(action, &QAction::triggered, this, [this, compAction]() {
                    // 1. 底层挂载组件
                    compAction.addComponent(m_Registry, m_CurrentEntity);
                    
                    // 2. 强制刷新 Inspector UI 
                    // (故意将 m_CurrentEntity 置空，绕过 BindEntity 开头的拦截)
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

        // 在按钮的正下方弹出菜单
        menu.exec(button->mapToGlobal(QPoint(0, button->height())));
    }

} // namespace Lizeral
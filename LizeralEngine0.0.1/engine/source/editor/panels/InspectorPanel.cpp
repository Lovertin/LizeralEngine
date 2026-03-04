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
        tryDrawComponent((ModelComponent*)nullptr,     "ModelComponent");

        // 在最底部加个弹簧，把所有组件往上顶
        m_MainLayout->addStretch(); 
    }
} // namespace Lizeral
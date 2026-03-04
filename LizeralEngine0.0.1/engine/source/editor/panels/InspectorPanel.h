#pragma once
#include <QWidget>
#include <QVBoxLayout>
#include <QGroupBox>
#include <QLabel>
#include "runtime/core/meta/reflection/reflection.h"
#include "runtime/function/ecs/entity.h"
#include "runtime/function/ecs/registry.h"

namespace Lizeral {
    class InspectorPanel : public QWidget {
        Q_OBJECT
    public:
        explicit InspectorPanel(QWidget* parent = nullptr);

        // 注入当前的 ECS 注册表（通常在编辑器启动时调用一次）
        void SetRegistry(Lizeral::Registry* registry) { m_Registry = registry; }

        void ClearPanel();

    public slots: // 【关键】：声明为 Qt 的槽函数，方便绑定全局事件
        // 当实体被选中时触发
        void BindEntity(Lizeral::Entity entity);

    private:
        QWidget* m_ContentWidget = nullptr; // 新增：用于包裹所有内容的容器
        QVBoxLayout* m_MainLayout;
        Lizeral::Registry* m_Registry { nullptr };
        Lizeral::Entity m_CurrentEntity { Lizeral::null_entity };
    };
}
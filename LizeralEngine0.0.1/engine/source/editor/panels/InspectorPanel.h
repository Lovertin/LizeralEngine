#pragma once 
// InspectorPanel.h
#pragma once
#include <QWidget>
#include <QVBoxLayout>

#include "editor/factory/ComponentUIFactory.h"
#include "runtime/function/ecs/entity.h"

namespace Lizeral{

    class InspectorPanel : public QWidget {
        Q_OBJECT
    public:
        explicit InspectorPanel(QWidget* parent = nullptr);
        
        // 核心接口：绑定当前选中的实体
        void BindEntity(Entity* entity);

    private:
        void ClearPanel(); // 清空当前面板上的所有组件 UI

        QVBoxLayout* m_MainLayout;
        Entity* m_CurrentEntity;
    };

}
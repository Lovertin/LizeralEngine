#include "InspectorPanel.h"
// #include "runtime/function/ecs/entity.h"

namespace Lizeral{

    InspectorPanel::InspectorPanel(QWidget* parent) : QWidget(parent), m_CurrentEntity(nullptr) {
        m_MainLayout = new QVBoxLayout(this);
        m_MainLayout->setAlignment(Qt::AlignTop);
        m_MainLayout->setSpacing(10);
    }

    void InspectorPanel::BindEntity(Entity* entity) {
        if (m_CurrentEntity == entity) return;
        
        m_CurrentEntity = entity;
        ClearPanel(); 

        if (!m_CurrentEntity) return;


        QGroupBox* testBox = new QGroupBox("Entity Header", this);
        QVBoxLayout* testLayout = new QVBoxLayout(testBox);
        testLayout->addWidget(new QLabel("Entity ID: 42", testBox));
        testLayout->addWidget(new QLabel("Name: Test GameObject", testBox));
        m_MainLayout->addWidget(testBox);
        
        QWidget* dummyComp = ComponentUIFactory::CreateWidgetFor(nullptr, this); // 触发兜底UI
        m_MainLayout->addWidget(dummyComp);
        // ----------------------------------------
    }

    // 必须加上这个空实现，否则会报错
    void InspectorPanel::ClearPanel() {
        // 如果 Layout 为空，直接返回
        if (!m_MainLayout) return;

        QLayoutItem* item;
        // 循环取出 Layout 中的第一个元素，直到取空
        while ((item = m_MainLayout->takeAt(0)) != nullptr) {
            // 如果这个元素包含 Widget
            if (QWidget* widget = item->widget()) {
                // 安全删除 Widget
                widget->deleteLater(); 
            }
            // 如果这个元素包含子 Layout，你也需要递归清理（当前简单起见先不写递归）
            
            // 删除 LayoutItem 本身
            delete item;
        }
    }

}
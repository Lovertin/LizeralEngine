#include "SceneOutlinerPanel.h"
#include <QVariant>
#include <QTimer>

// 请根据你的实际路径修改以下两个头文件
#include "editor/selection/EditorSelection.h" 
#include "runtime/function/ecs/components/Name/NameComponent.h" 
#include "runtime/function/ecs/components/EditorOnly/EditorOnlyComponent.h"\

#include <iostream>

namespace Lizeral {

    SceneOutlinerPanel::SceneOutlinerPanel(QWidget* parent) : QWidget(parent) {
        m_MainLayout = new QVBoxLayout(this);
        m_MainLayout->setContentsMargins(0, 0, 0, 0); // 去掉边距，让树形控件填满 Dock

        m_TreeWidget = new QTreeWidget(this);
        m_TreeWidget->setHeaderLabel("Hierarchy");
        m_TreeWidget->setContextMenuPolicy(Qt::CustomContextMenu); // 允许自定义右键菜单
        
        m_MainLayout->addWidget(m_TreeWidget);

        // 绑定事件
        connect(m_TreeWidget, &QTreeWidget::itemClicked, this, &SceneOutlinerPanel::OnItemClicked);
        connect(m_TreeWidget, &QTreeWidget::customContextMenuRequested, this, &SceneOutlinerPanel::ShowContextMenu);
    }

    void SceneOutlinerPanel::SetRegistry(Lizeral::Registry* registry) {
        m_Registry = registry;
        Refresh(); // 注入成功后立刻刷新一次
    }

    void SceneOutlinerPanel::Refresh() {
        // std::cout<<"start cleaning TreeWidget"<<std::endl;
        m_TreeWidget->clear();
        // std::cout<<"m_TreeWidget is cleared"<<std::endl;
        if (!m_Registry) return;

        // 遍历所有拥有 NameComponent 的实体
        auto view = m_Registry->view<NameComponent>();
        for (auto entity : view) {
            // 这里你可以加上黑名单检查（比如检查有没有 EditorOnlyComponent）
            if(m_Registry->has<EditorOnlyComponent>(entity)){
                continue;
            }
            
            auto& nameComp = m_Registry->get<NameComponent>(entity);
            
            QTreeWidgetItem* item = new QTreeWidgetItem(m_TreeWidget);
            item->setText(0, QString::fromStdString(nameComp.getName()));
            
            // 【核心魔法】：把 Entity ID 藏在 Item 的用户数据区 (UserRole)
            item->setData(0, Qt::UserRole, QVariant::fromValue(static_cast<uint32_t>(entity)));
        }
    }

    void SceneOutlinerPanel::OnItemClicked(QTreeWidgetItem* item, int column) {
        if (!item) return;

        // 取出藏在里面的 Entity ID
        uint32_t entityId = item->data(0, Qt::UserRole).toUInt();
        Lizeral::Entity entity = static_cast<Lizeral::Entity>(entityId);

        // 通知全局系统：我选中了这个实体！
        EditorSelection::Get().SelectEntity(entity);
    }

    void SceneOutlinerPanel::ShowContextMenu(const QPoint& pos) {
        if (!m_Registry) return;

        QMenu contextMenu(this);
        QAction* createEntityAction = contextMenu.addAction("Create Empty Entity");

        // 在鼠标点击的全局屏幕位置弹出菜单，并阻塞等待用户选择
        QAction* selectedAction = contextMenu.exec(m_TreeWidget->mapToGlobal(pos));

        if (selectedAction == createEntityAction) {
            // 使用 Qt::QueuedConnection 确保完全推迟到下一帧的主线程事件循环中执行
            QMetaObject::invokeMethod(this, [this]() {
                // std::cout << "Creating new entity in main thread..." << std::endl;
                Lizeral::Entity newEntity = m_Registry->create();
                m_Registry->emplace<NameComponent>(newEntity, "New Entity");
                
                // std::cout << "Starting Refresh..." << std::endl;
                Refresh(); // 此时再清理树节点
                
                EditorSelection::Get().SelectEntity(newEntity);
            }, Qt::QueuedConnection);
        }
    }

} // namespace Lizeral
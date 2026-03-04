#include "SceneOutlinerPanel.h"
#include <QVariant>
#include <QTimer>

// 请根据你的实际路径修改以下两个头文件
#include "editor/selection/EditorSelection.h" 
#include "runtime/function/ecs/components/componentAll.h"

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
        contextMenu.setStyleSheet("QMenu { background-color: #333333; color: white; border: 1px solid #555555; } QMenu::item:selected { background-color: #4CAF50; }");

        // ==========================================
        // 第一步：先检测鼠标点在了哪里，并组装所有菜单项
        // ==========================================
        QAction* deleteEntityAction = nullptr;
        QTreeWidgetItem* clickedItem = m_TreeWidget->itemAt(pos);
        
        // 1. 如果点中了一个实体，增加“删除”选项
        if (clickedItem) {
            deleteEntityAction = contextMenu.addAction("Delete Entity");
            // 这里我们用样式表覆盖一下，让删除按钮变红
            contextMenu.setStyleSheet(contextMenu.styleSheet() + " QAction#DeleteAction { color: #ff6666; font-weight: bold; }");
            deleteEntityAction->setObjectName("DeleteAction");
            contextMenu.addSeparator(); // 加一条分割线
        }

        // 2. 组装“创建”预设实体的菜单项
        QAction* createEmptyAction = contextMenu.addAction("Create Empty Entity");
        QMenu* object3DMenu = contextMenu.addMenu("3D Object");
        QAction* createCubeAction = object3DMenu->addAction("Cube");
        QMenu* lightMenu = contextMenu.addMenu("Light");
        QAction* createDirLightAction = lightMenu->addAction("Directional Light");
        QAction* createCameraAction = contextMenu.addAction("Camera");

        // ==========================================
        // 第二步：将菜单显示在屏幕上，并阻塞等待用户点击
        // ==========================================
        QAction* selectedAction = contextMenu.exec(m_TreeWidget->mapToGlobal(pos));

        // ==========================================
        // 第三步：根据用户的点击，执行对应的底层逻辑
        // ==========================================
        if (!selectedAction) return; // 用户什么都没点，点击了旁边取消了菜单

        if (selectedAction == deleteEntityAction && clickedItem) {
            // --- 处理实体删除 ---
            uint32_t entityId = clickedItem->data(0, Qt::UserRole).toUInt();
            Lizeral::Entity targetEntity = static_cast<Lizeral::Entity>(entityId);

            QMessageBox::StandardButton reply = QMessageBox::question(this, "Delete Entity", 
                "Delete '" + clickedItem->text(0) + "' forever?", QMessageBox::Yes | QMessageBox::No);
            
            if (reply == QMessageBox::Yes) {
                // 如果删除的是当前正选中的实体，要清空右侧 Inspector 面板
                if (EditorSelection::Get().GetSelected() == targetEntity) {
                    EditorSelection::Get().SelectEntity(Lizeral::null_entity);
                }
                
                // 1. 底层彻底销毁实体（调用你刚写的 destroy）
                m_Registry->destroy(targetEntity);
                
                // 2. 刷新左侧大纲树
                Refresh();
            }
        } 
        else {
            // --- 处理实体创建 ---
            QMetaObject::invokeMethod(this, [this, selectedAction, createEmptyAction, createCubeAction, createDirLightAction, createCameraAction]() {
                
                Lizeral::Entity newEntity = m_Registry->create();
                m_Registry->emplace<TransformComponent>(newEntity); // 都有 Transform

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

                Refresh(); // 刷新左侧树
                EditorSelection::Get().SelectEntity(newEntity); // 通知右侧 Inspector 选中它

            }, Qt::QueuedConnection);
        }
    }

} // namespace Lizeral
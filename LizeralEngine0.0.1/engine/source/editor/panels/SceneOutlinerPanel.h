#pragma once
#include <QWidget>
#include <QTreeWidget>
#include <QVBoxLayout>
#include <QMenu>
#include <QMessageBox>
#include <QAction>

// 包含你的 ECS 核心
#include "runtime/function/ecs/entity.h"
#include "runtime/function/ecs/registry.h"

namespace Lizeral {

    class SceneOutlinerPanel : public QWidget {
        Q_OBJECT
    public:
        explicit SceneOutlinerPanel(QWidget* parent = nullptr);

        // 注入全局 Registry
        void SetRegistry(Lizeral::Registry* registry);

        // 刷新列表（当新建、删除实体时调用）
        void Refresh();

    private slots:
        // Qt 槽函数：处理点击和右键菜单
        void OnItemClicked(QTreeWidgetItem* item, int column);
        void ShowContextMenu(const QPoint& pos);

    private:
        QVBoxLayout* m_MainLayout;
        QTreeWidget* m_TreeWidget;
        Lizeral::Registry* m_Registry { nullptr };
    };

} // namespace Lizeral
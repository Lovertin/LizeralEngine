#pragma once
#include <QWidget>
#include <QTreeWidget>
#include <QVBoxLayout>
#include <QMenu>
#include <QMessageBox>
#include <QAction>

#include "runtime/function/ecs/entity.h"
#include "runtime/function/ecs/registry.h"

namespace Lizeral {

    class SceneOutlinerPanel : public QWidget {
        Q_OBJECT
    public:
        explicit SceneOutlinerPanel(QWidget* parent = nullptr);

        void SetRegistry(Lizeral::Registry* registry);

        void Refresh();

    private slots:
        void OnItemClicked(QTreeWidgetItem* item, int column);
        void ShowContextMenu(const QPoint& pos);

    private:
        QVBoxLayout* m_MainLayout;
        QTreeWidget* m_TreeWidget;
        Lizeral::Registry* m_Registry { nullptr };
    };

} // namespace Lizeral
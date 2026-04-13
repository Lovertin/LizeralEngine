#pragma once
#include <QWidget>
#include <QVBoxLayout>
#include <QGroupBox>
#include <QLabel>
#include <QPushButton>
#include <QMenu>
#include <QMessageBox>

#include <functional>
#include <vector>
#include <string>

#include "runtime/core/meta/reflection/reflection.h"
#include "runtime/function/ecs/entity.h"
#include "runtime/function/ecs/registry.h"

namespace Lizeral {

    struct ComponentAction {
        std::string name;
        std::function<bool(Lizeral::Registry*, Lizeral::Entity)> hasComponent;
        std::function<void(Lizeral::Registry*, Lizeral::Entity)> addComponent;
    };

    class InspectorPanel : public QWidget {
        Q_OBJECT
    public:
        explicit InspectorPanel(QWidget* parent = nullptr);

        void SetRegistry(Lizeral::Registry* registry) { m_Registry = registry; }

        void ClearPanel();

    public slots: 
        void BindEntity(Lizeral::Entity entity);

    private:

        template<typename T>
        void RegisterComponentType(const std::string& name) {
            ComponentAction action;
            action.name = name;
            action.hasComponent = [](Lizeral::Registry* r, Lizeral::Entity e) { return r->has<T>(e); };
            action.addComponent = [](Lizeral::Registry* r, Lizeral::Entity e) { r->emplace<T>(e); };
            m_AvailableComponents.push_back(action);
        }

        void InitComponentRegistry();

        void ShowAddComponentMenu(QPushButton* button);

    private:    

        QWidget* m_ContentWidget = nullptr;
        QVBoxLayout* m_MainLayout;
        Registry* m_Registry { nullptr };
        Entity m_CurrentEntity { Lizeral::null_entity };
        int32_t m_CurrentMeshAssetIndex { -1 };

        std::vector<ComponentAction> m_AvailableComponents;
    };
}

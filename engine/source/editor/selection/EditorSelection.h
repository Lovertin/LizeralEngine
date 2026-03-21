// EditorSelection.h
#pragma once
#include <QObject>
#include "runtime/function/ecs/entity.h"

namespace Lizeral {

    class EditorSelection : public QObject {
        Q_OBJECT
    public:
        // singleton class
        static EditorSelection& Get() {
            static EditorSelection* instance = new EditorSelection();
            return *instance;
        }

        // Disable copy constructor and move constructor
        EditorSelection(const EditorSelection&) = delete;
        EditorSelection& operator=(const EditorSelection&) = delete;
        EditorSelection(EditorSelection&&) = delete;
        EditorSelection& operator=(EditorSelection&&) = delete;

        void SelectEntity(Lizeral::Entity entity);

        Lizeral::Entity GetSelected() const { return m_SelectedEntity; }

    signals:
        void OnEntitySelected(Lizeral::Entity entity);

        void OnEntityDataModified(Lizeral::Entity entity);

    private:

        explicit EditorSelection(QObject* parent = nullptr) : QObject(parent) {}
        
        ~EditorSelection() override = default;

    private:
        Lizeral::Entity m_SelectedEntity{ Lizeral::null_entity };
    };
    
} // namespace Lizeral
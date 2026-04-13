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
        void SelectSubMesh(Lizeral::Entity entity, uint32_t mesh_asset_index);
        void ClearSubMeshSelection();

        Lizeral::Entity GetSelected() const { return m_SelectedEntity; }
        int32_t GetSelectedSubMeshIndex() const { return m_SelectedMeshAssetIndex; }
        bool HasSelectedSubMesh() const { return m_SelectedMeshAssetIndex >= 0; }
        bool IsSubMeshSelectedForEntity(Lizeral::Entity entity) const {
            return m_SelectedEntity == entity && m_SelectedMeshAssetIndex >= 0;
        }

    signals:
        void OnEntitySelected(Lizeral::Entity entity);
        void OnSubMeshSelectionChanged(Lizeral::Entity entity, int32_t meshAssetIndex);

        void OnEntityDataModified(Lizeral::Entity entity);

    private:

        explicit EditorSelection(QObject* parent = nullptr) : QObject(parent) {}
        
        ~EditorSelection() override = default;

    private:
        Lizeral::Entity m_SelectedEntity{ Lizeral::null_entity };
        int32_t m_SelectedMeshAssetIndex { -1 };
    };
    
} // namespace Lizeral

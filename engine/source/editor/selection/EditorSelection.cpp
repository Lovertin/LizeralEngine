#include "EditorSelection.h"

namespace Lizeral{

    void EditorSelection::SelectEntity(Lizeral::Entity entity)
    {
            const bool entityChanged = (m_SelectedEntity != entity);
            const bool meshChanged = (m_SelectedMeshAssetIndex != -1);

            if (!entityChanged && !meshChanged) {
                return;
            }

            m_SelectedEntity = entity;
            m_SelectedMeshAssetIndex = -1;

            if (entityChanged) {
                emit OnEntitySelected(entity);
            }
            emit OnSubMeshSelectionChanged(entity, -1);
    }

    void EditorSelection::SelectSubMesh(Lizeral::Entity entity, uint32_t mesh_asset_index) {
            const int32_t meshIndex = static_cast<int32_t>(mesh_asset_index);
            const bool entityChanged = (m_SelectedEntity != entity);
            const bool meshChanged = (m_SelectedMeshAssetIndex != meshIndex);

            if (!entityChanged && !meshChanged) {
                return;
            }

            m_SelectedEntity = entity;
            m_SelectedMeshAssetIndex = meshIndex;

            if (entityChanged) {
                emit OnEntitySelected(entity);
            }
            emit OnSubMeshSelectionChanged(entity, meshIndex);
    }

    void EditorSelection::ClearSubMeshSelection() {
            if (m_SelectedMeshAssetIndex == -1) {
                return;
            }

            m_SelectedMeshAssetIndex = -1;
            emit OnSubMeshSelectionChanged(m_SelectedEntity, -1);
    }
}


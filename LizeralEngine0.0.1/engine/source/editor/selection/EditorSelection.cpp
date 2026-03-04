#include "EditorSelection.h"

namespace Lizeral{

    void EditorSelection::SelectEntity(Lizeral::Entity entity)
    {
            if (m_SelectedEntity != entity) {
                m_SelectedEntity = entity;
                emit OnEntitySelected(entity); // 广播给全宇宙！
            }
    }
}


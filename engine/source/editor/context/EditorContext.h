#pragma once
#include "editor/command/CommandManager.h"
#include "runtime/function/ecs/registry.h"
#include "runtime/function/ecs/entity.h"
#include "editor/selection/EditorSelection.h"
#include "runtime/function/ecs/components/componentAll.h"
#include <string>

namespace Lizeral {
    enum class EditorMode{
        Edit,
        Play
    };

    class EditorContext {
    public:
        static EditorContext& Get() {
            static EditorContext instance;
            return instance;
        }

        CommandManager* GetCommandManager() { return &m_CommandManager; }
        
        void SetRegistry(Registry* reg) { m_Registry = reg; }
        Registry* GetRegistry() { return m_Registry; }

        void* GetComponentByName(Entity entity, const std::string& compName);

        EditorMode getMode() const { return m_mode ;} 

        void setMode(const EditorMode mode) {m_mode = mode ;}

    private:
        EditorContext() = default;
        CommandManager m_CommandManager;
        Registry* m_Registry = nullptr;
        EditorMode m_mode;

    };

} // namespace Lizeral
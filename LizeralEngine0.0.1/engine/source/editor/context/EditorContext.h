#pragma once
#include "editor/command/CommandManager.h"
#include "runtime/function/ecs/registry.h"
#include "runtime/function/ecs/entity.h"
#include "editor/selection/EditorSelection.h"
#include "runtime/function/ecs/components/componentAll.h"
#include <string>

namespace Lizeral {

    class EditorContext {
    public:
        static EditorContext& Get() {
            static EditorContext instance;
            return instance;
        }

        CommandManager* GetCommandManager() { return &m_CommandManager; }
        
        void SetRegistry(Registry* reg) { m_Registry = reg; }
        Registry* GetRegistry() { return m_Registry; }

        // 将反射系统的字符串，翻译为 ECS 的真实内存指针！
        void* GetComponentByName(Entity entity, const std::string& compName);

    private:
        EditorContext() = default;
        CommandManager m_CommandManager;
        Registry* m_Registry = nullptr;
    };

} // namespace Lizeral
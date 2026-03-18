#pragma once
#include "runtime/function/ecs/registry.h"
#include "runtime/core/meta/serializer/serializer.h"
#include <string>

namespace Lizeral {

    class SceneSerializer {
    public:
        // 1. 内存级：将当前 Registry 序列化为 JSON 对象 (用于 Play/Stop 快照)
        static PJson Serialize(Registry* registry);

        // 2. 内存级：清空当前 Registry，并从 JSON 对象恢复场景
        static void Deserialize(const PJson& json, Registry* registry);

        // 3. 磁盘级：保存场景到本地文件 (例如 File -> Save)
        static void SaveToFile(Registry* registry, const std::string& filepath);

        // 4. 磁盘级：从本地文件加载场景 (例如 File -> Load)
        static void LoadFromFile(const std::string& filepath, Registry* registry);
    };

} // namespace Lizeral
#pragma once
#include "runtime/function/ecs/registry.h"
#include "runtime/core/meta/serializer/serializer.h"
#include <sstream>
#include <string>

namespace Lizeral {

    class SceneSerializer {
    public:
        static PJson Serialize(Registry* registry);

        static void Deserialize(const PJson& json, Registry* registry);

        static void SaveToFile(Registry* registry, const std::string& filepath);

        static void LoadFromFile(const std::string& filepath, Registry* registry);
    };

} // namespace Lizeral
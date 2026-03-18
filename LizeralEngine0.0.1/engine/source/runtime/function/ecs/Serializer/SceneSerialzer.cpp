#include "SceneSerializer.h"
#include "runtime/function/ecs/components/componentAll.h"
#include "_generated/serializer/all_serializer.ipp"
#include <iostream>
#include <fstream>

namespace Lizeral {

    // --- 两个极其好用的内部模板辅助函数 ---

    // 辅助函数：如果实体有这个组件，就把它转成 JSON
    template<typename T>
    static void TrySerializeComponent(Registry* registry, Entity entity, PJson::object& entityJson, const std::string& compName) {
        if (registry->has<T>(entity)) {
            // 调用你自动生成的 PSerializer::write
            entityJson[compName] = PSerializer::write(registry->get<T>(entity));
        }
    }

    // 辅助函数：如果 JSON 里有这个组件的数据，就给实体挂载并反序列化
    template<typename T>
    static void TryDeserializeComponent(Registry* registry, Entity entity, const PJson& entityJson, const std::string& compName) {
        if (!entityJson[compName].is_null()) {
            auto& comp = registry->emplace<T>(entity);
            // 调用你自动生成的 PSerializer::read
            PSerializer::read(entityJson[compName], comp);
        }
    }

    // ==========================================
    // 核心实现
    // ==========================================

    PJson SceneSerializer::Serialize(Registry* registry) {
        PJson::array sceneArray;

        auto view = registry->view<NameComponent>();
        for (auto entity : view) {
            PJson::object entityJson;
            
            // 记录实体的原始 ID（可选，主要用于调试）
            // entityJson["EntityID"] = PSerializer::write(static_cast<uint32_t>(entity));

            // 依次检查并序列化所有已知的组件
            TrySerializeComponent<NameComponent>(registry, entity, entityJson, "NameComponent");
            TrySerializeComponent<TransformComponent>(registry, entity, entityJson, "TransformComponent");
            TrySerializeComponent<RigidBodyComponent>(registry, entity, entityJson, "RigidBodyComponent");
            TrySerializeComponent<ColliderComponent>(registry, entity, entityJson, "ColliderComponent");
            TrySerializeComponent<CameraComponent>(registry, entity, entityJson, "CameraComponent");
            TrySerializeComponent<CameraControlComponent>(registry, entity, entityJson, "CameraControlComponent");
            TrySerializeComponent<DirectionLightComponent>(registry, entity, entityJson, "DirectionLightComponent");
            TrySerializeComponent<VulkanModelComponent>(registry, entity, entityJson, "VulkanModelComponent");

            // 将这个实体的 JSON 对象加入场景数组
            sceneArray.push_back(PJson(entityJson));
        }

        PJson::object root;
        root["Entities"] = sceneArray;
        return PJson(root);
    }

    void SceneSerializer::Deserialize(const PJson& json, Registry* registry) {
        if (json.is_null() || json["Entities"].is_null()) {
            std::cout << "[SceneSerializer] Failed to load: JSON is empty or missing 'Entities' array." << std::endl;
            return;
        }

        // 1. 彻底清空当前世界！
        registry->clear();

        // 2. 遍历 JSON 数组，重建世界
        const auto& sceneArray = json["Entities"].array_items();
        for (const auto& entityItem : sceneArray) {
            // 创建一个全新的干净实体
            Entity newEntity = registry->create();

            // 依次尝试反序列化各个组件
            TryDeserializeComponent<NameComponent>(registry, newEntity, entityItem, "NameComponent");
            TryDeserializeComponent<TransformComponent>(registry, newEntity, entityItem, "TransformComponent");
            TryDeserializeComponent<RigidBodyComponent>(registry, newEntity, entityItem, "RigidBodyComponent");
            TryDeserializeComponent<ColliderComponent>(registry, newEntity, entityItem, "ColliderComponent");
            TryDeserializeComponent<CameraComponent>(registry, newEntity, entityItem, "CameraComponent");
            TryDeserializeComponent<CameraControlComponent>(registry, newEntity, entityItem, "CameraControlComponent");
            TryDeserializeComponent<DirectionLightComponent>(registry, newEntity, entityItem, "DirectionLightComponent");
            TryDeserializeComponent<VulkanModelComponent>(registry, newEntity, entityItem, "VulkanModelComponent");
        }
    }

    // --- 磁盘操作 ---
    void SceneSerializer::SaveToFile(Registry* registry, const std::string& filepath) {
        // 你可以根据你的 PJson 库提供的 dump() 序列化为字符串
        // PJson json = Serialize(registry);
        // std::ofstream out(filepath);
        // out << json.dump(); 
        // out.close();
    }

    void SceneSerializer::LoadFromFile(const std::string& filepath, Registry* registry) {
        // 读取文本文件，解析为 PJson，然后调用 Deserialize()
    }

} // namespace Lizeral
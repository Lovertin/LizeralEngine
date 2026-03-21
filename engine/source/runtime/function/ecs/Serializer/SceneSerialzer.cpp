#include "SceneSerializer.h"
#include "runtime/function/ecs/components/componentAll.h"
#include "_generated/serializer/all_serializer.ipp"
#include <iostream>
#include <fstream>

namespace Lizeral {

    // trasform to JSON
    template<typename T>
    static void TrySerializeComponent(Registry* registry, Entity entity, PJson::object& entityJson, const std::string& compName) {
        if (registry->has<T>(entity)) {
            // PSerializer::write
            entityJson[compName] = PSerializer::write(registry->get<T>(entity));
        }
    }

    template<typename T>
    static void TryDeserializeComponent(Registry* registry, Entity entity, const PJson& entityJson, const std::string& compName) {
        if (!entityJson[compName].is_null()) {
            auto& comp = registry->emplace<T>(entity);
            // PSerializer::read
            PSerializer::read(entityJson[compName], comp);
        }
    }

    PJson SceneSerializer::Serialize(Registry* registry) {
        PJson::array sceneArray;

        auto view = registry->view<NameComponent>();
        for (auto entity : view) {
            PJson::object entityJson;

            TrySerializeComponent<NameComponent>(registry, entity, entityJson, "NameComponent");
            TrySerializeComponent<TransformComponent>(registry, entity, entityJson, "TransformComponent");
            TrySerializeComponent<RigidBodyComponent>(registry, entity, entityJson, "RigidBodyComponent");
            TrySerializeComponent<ColliderComponent>(registry, entity, entityJson, "ColliderComponent");
            TrySerializeComponent<CameraComponent>(registry, entity, entityJson, "CameraComponent");
            TrySerializeComponent<CameraControlComponent>(registry, entity, entityJson, "CameraControlComponent");
            TrySerializeComponent<DirectionLightComponent>(registry, entity, entityJson, "DirectionLightComponent");
            TrySerializeComponent<PointLightComponent>(registry, entity, entityJson, "PointLightComponent");
            TrySerializeComponent<VulkanModelComponent>(registry, entity, entityJson, "VulkanModelComponent");

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

        registry->clear();

        // rebuild world
        const auto& sceneArray = json["Entities"].array_items();
        for (const auto& entityItem : sceneArray) {

            Entity newEntity = registry->create();

            TryDeserializeComponent<NameComponent>(registry, newEntity, entityItem, "NameComponent");
            TryDeserializeComponent<TransformComponent>(registry, newEntity, entityItem, "TransformComponent");
            TryDeserializeComponent<RigidBodyComponent>(registry, newEntity, entityItem, "RigidBodyComponent");
            TryDeserializeComponent<ColliderComponent>(registry, newEntity, entityItem, "ColliderComponent");
            TryDeserializeComponent<CameraComponent>(registry, newEntity, entityItem, "CameraComponent");
            TryDeserializeComponent<CameraControlComponent>(registry, newEntity, entityItem, "CameraControlComponent");
            TryDeserializeComponent<DirectionLightComponent>(registry, newEntity, entityItem, "DirectionLightComponent");
            TryDeserializeComponent<PointLightComponent>(registry, newEntity, entityItem, "PointLightComponent");
            TryDeserializeComponent<VulkanModelComponent>(registry, newEntity, entityItem, "VulkanModelComponent");
        }
    }

    // --- 磁盘操作 ---
    void SceneSerializer::SaveToFile(Registry* registry, const std::string& filepath) {
        // 1. 获取场景的 JSON 对象
        PJson json = Serialize(registry);
        
        // 2. 使用 json11 的 dump() 转为格式化字符串
        std::string jsonString = json.dump(); 

        // 3. 写入文件
        std::ofstream out(filepath);
        if (out.is_open()) {
            out << jsonString;
            out.close();
            std::cout << "[SceneSerializer] Successfully saved scene to: " << filepath << std::endl;
        } else {
            std::cerr << "[SceneSerializer] Error: Could not open file for writing: " << filepath << std::endl;
        }
    }

    void SceneSerializer::LoadFromFile(const std::string& filepath, Registry* registry) {
        // 1. 读取文件内容
        std::ifstream in(filepath);
        if (!in.is_open()) {
            std::cerr << "[SceneSerializer] Error: Could not open file for reading: " << filepath << std::endl;
            return;
        }

        std::stringstream buffer;
        buffer << in.rdbuf();
        std::string jsonString = buffer.str();
        in.close();

        // 2. 解析 JSON 字符串
        std::string err;
        PJson json = PJson::parse(jsonString, err);
        
        if (!err.empty()) {
            std::cerr << "[SceneSerializer] JSON Parsing Error: " << err << std::endl;
            return;
        }

        // 3. 反序列化，重建 ECS 世界！
        Deserialize(json, registry);
        std::cout << "[SceneSerializer] Successfully loaded scene from: " << filepath << std::endl;
    }

} // namespace Lizeral
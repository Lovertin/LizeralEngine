#include "EditorContext.h"

namespace Lizeral{

    void* EditorContext::GetComponentByName(Entity entity, const std::string& compName) {
            if (!m_Registry || entity == null_entity) return nullptr;

            // 这里做一个简单的路由（后期可以用宏或模板自动化）
            if (compName == "TransformComponent" && m_Registry->has<TransformComponent>(entity)) 
                return &m_Registry->get<TransformComponent>(entity);
            if (compName == "RigidBodyComponent" && m_Registry->has<RigidBodyComponent>(entity)) 
                return &m_Registry->get<RigidBodyComponent>(entity);
            if (compName == "ColliderComponent" && m_Registry->has<ColliderComponent>(entity)) 
                return &m_Registry->get<ColliderComponent>(entity);
            if (compName == "DirectionLightComponent" && m_Registry->has<DirectionLightComponent>(entity)) 
                return &m_Registry->get<DirectionLightComponent>(entity);
            if (compName == "CameraComponent" && m_Registry->has<CameraComponent>(entity)) 
                return &m_Registry->get<CameraComponent>(entity);
            if (compName == "VulkanModelComponent" && m_Registry->has<VulkanModelComponent>(entity))
                return &m_Registry->get<VulkanModelComponent>(entity);
            if (compName == "CameraControlComponent" && m_Registry->has<CameraControlComponent>(entity))
                return &m_Registry->get<CameraControlComponent>(entity);

            return nullptr;
        }

}
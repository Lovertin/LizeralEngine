#pragma once 
#include "runtime/core/meta/reflection/reflection.h"
#include "runtime/function/ecs/components/component.h"
#include <string>

namespace Lizeral{
        
    REFLECTION_TYPE(VulkanModelComponent)
    CLASS(VulkanModelComponent : public Component, WhiteListFields){
        REFLECTION_BODY(VulkanModelComponent)
    public:
        std::string getVulkanModelPath() const { return m_ModelPath; }
        
        void setVulkanModelPath(const std::string& glb_path) { 
            if (m_ModelPath != glb_path) {
                m_ModelPath = glb_path;
                m_IsLoaded = false; 
            }
        } 

        bool IsLoaded() const { return m_IsLoaded; }
        void SetLoaded(bool state) { m_IsLoaded = state; }

    private:
        META(Enable, UI:Address)
        std::string m_ModelPath;

        bool m_IsLoaded = false; 
    };

}
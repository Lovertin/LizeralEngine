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
            // 只有当路径真正发生改变时，才重置状态
            if (m_ModelPath != glb_path) {
                m_ModelPath = glb_path;
                m_IsLoaded = false; //  脏标记：告诉系统我换模型了，需要重新加载！
            }
        } 

        // 提供给 Render System 检查的状态接口
        bool IsLoaded() const { return m_IsLoaded; }
        void SetLoaded(bool state) { m_IsLoaded = state; }

    private:
        META(Enable, UI:Address)
        std::string m_ModelPath;

        //运行时状态 (不加 META，不需要序列化，也不需要显示在 UI 上)
        bool m_IsLoaded = false; 
    };

}
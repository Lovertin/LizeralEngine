#pragma once

#include <memory>
#include <string>
#include "runtime/function/res_type/shader/shader.h"
#include "runtime/core/meta/reflection/reflection.h" // 【新增】：引入反射宏

namespace Lizeral {

    REFLECTION_TYPE(Material)
    CLASS(Material, WhiteListFields) {
        REFLECTION_BODY(Material)
    protected:
        // ==========================================
        // 【运行时数据】：不反射，只在内存中起作用
        // ==========================================
        std::shared_ptr<Shader> m_Shader;

    public:
        // ==========================================
        // 【反射数据】：供 UI 显示和场景保存的路径
        // ==========================================
        META(Enable)
        std::string m_VertShaderPath;
        
        META(Enable)
        std::string m_FragShaderPath;

        Material() = default;
        Material(std::shared_ptr<Shader> shader) : m_Shader(shader) {}
        
        virtual ~Material() = default;

        // 核心接口
        virtual void BindAndApply(){}

        // 【新增】：根据反射的路径，真正去加载 Shader
        // 可以放在 cpp 里实现，需要调用 ResourceManager
        void LoadShader(); 

        // Getter / Setter
        void SetShader(std::shared_ptr<Shader> shader) { m_Shader = shader; }
        std::shared_ptr<Shader> GetShader() const { return m_Shader; }
    };

}
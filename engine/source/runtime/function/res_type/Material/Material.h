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

        std::shared_ptr<Shader> m_Shader;

    public:

        META(Enable)
        std::string m_VertShaderPath;
        
        META(Enable)
        std::string m_FragShaderPath;

        Material() = default;
        Material(std::shared_ptr<Shader> shader) : m_Shader(shader) {}
        
        virtual ~Material() = default;

        virtual void BindAndApply(){}

        void LoadShader(); 

        // Getter / Setter
        void SetShader(std::shared_ptr<Shader> shader) { m_Shader = shader; }
        std::shared_ptr<Shader> GetShader() const { return m_Shader; }
    };

}
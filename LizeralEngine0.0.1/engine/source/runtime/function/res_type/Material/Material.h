// Material.h
#pragma once

#include <memory>
#include "runtime/function/res_type/shader/shader.h"

namespace Lizeral {

    class Material {
    protected:
        // 每个材质实例都必须绑定一个具体的 Shader
        std::shared_ptr<Shader> m_Shader;

    public:
        Material() = default;
        Material(std::shared_ptr<Shader> shader) : m_Shader(shader) {}
        
        // 虚析构函数，保证派生类被正确清理
        virtual ~Material() = default;

        // ==========================================
        // 核心接口：子类必须实现该方法。
        // 在这里执行 shader->Bind()，并把特有的 Uniform 送进 GPU。
        // ==========================================
        virtual void BindAndApply() = 0;

        // Getter / Setter
        void SetShader(std::shared_ptr<Shader> shader) { m_Shader = shader; }
        std::shared_ptr<Shader> GetShader() const { return m_Shader; }
    };

}
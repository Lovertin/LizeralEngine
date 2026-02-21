// PBRMaterial.h
#pragma once

#include "Material.h"
#include "runtime/core/math/vector3.h"
// 预留：#include "runtime/resource/texture/texture2d.h" 

namespace Lizeral {

    class PBRMaterial : public Material {
    public:
        // ------------------------------------------
        // PBR 核心参数 (基础版：纯数值驱动)
        // ------------------------------------------
        Vector3 m_Albedo { 1.0f, 1.0f, 1.0f }; // 反照率 (基础颜色)
        float   m_Metallic { 0.0f };           // 金属度：0是绝缘体(塑料/木头)，1是金属
        float   m_Roughness { 0.5f };          // 粗糙度：0是绝对光滑(镜子)，1是绝对粗糙(泥土)
        float   m_AO { 1.0f };                 // 环境光遮蔽 (Ambient Occlusion)

        // ------------------------------------------
        // PBR 核心参数 (进阶版：贴图驱动，体现架构前瞻性)
        // std::shared_ptr<Texture2D> m_AlbedoMap;
        // std::shared_ptr<Texture2D> m_MetallicMap;
        // std::shared_ptr<Texture2D> m_RoughnessMap;
        // ------------------------------------------

        PBRMaterial() = default;
        PBRMaterial(std::shared_ptr<Shader> shader) : Material(shader) {}

        // 实现基类的纯虚函数
        void BindAndApply() override {
            if (!m_Shader) return;

            // 1. 激活 PBR 专属的 Shader
            m_Shader->Bind();

            // 2. 上传当前材质独有的物理参数到 GPU
            m_Shader->SetUniformVector3f("u_Albedo", m_Albedo);
            m_Shader->SetUniform1f("u_Metallic", m_Metallic);
            m_Shader->SetUniform1f("u_Roughness", m_Roughness);
            m_Shader->SetUniform1f("u_AO", m_AO);

            // 3. 如果未来加载了贴图，在这里绑定纹理单元并上传
            // if (m_AlbedoMap) {
            //     m_AlbedoMap->Bind(0); 
            //     m_Shader->SetUniform1i("u_AlbedoMap", 0);
            // }
        }
    };

}
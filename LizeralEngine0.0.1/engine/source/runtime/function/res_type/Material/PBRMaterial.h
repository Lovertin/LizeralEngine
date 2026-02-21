#pragma once

#include "Material.h"
#include "runtime/core/math/vector3.h"
#include "runtime/function/res_type/texture/Texture2D.h" // 【新增】：引入刚刚写好的 Texture2D
#include "runtime/function/res_type/texture/TextureCube.h"

namespace Lizeral {

    class PBRMaterial : public Material {
    public:
        // PBR 数值参数
        Vector3 m_Albedo { 1.0f, 1.0f, 1.0f };
        float   m_Metallic { 0.0f };          
        float   m_Roughness { 0.5f };         
        float   m_AO { 1.0f };                

        // 【解封】：PBR 贴图指针
        std::shared_ptr<Texture2D> m_AlbedoMap;
        std::shared_ptr<TextureCube> m_IrradianceMap;
        std::shared_ptr<TextureCube> m_PrefilterMap;

        PBRMaterial() = default;
        PBRMaterial(std::shared_ptr<Shader> shader) : Material(shader) {}

        void BindAndApply() override {
            if (!m_Shader) return;

            m_Shader->Bind();

            // 上传数值
            m_Shader->SetUniformVector3f("u_Albedo", m_Albedo);
            m_Shader->SetUniform1f("u_Metallic", m_Metallic);
            m_Shader->SetUniform1f("u_Roughness", m_Roughness);
            m_Shader->SetUniform1f("u_AO", m_AO);

            // 【核心挂载逻辑】：判断是否有贴图，并正确绑定纹理槽位
            if (m_AlbedoMap) {
                m_AlbedoMap->Bind(0); // 绑定到显卡的 0 号纹理插槽
                m_Shader->SetUniform1i("u_AlbedoMap", 0); // 告诉 Shader 采样器去 0 号槽读图
                m_Shader->SetUniform1i("u_UseAlbedoMap", 1); // 开启 Shader 里的贴图开关 (true)
            } else {
                m_Shader->SetUniform1i("u_UseAlbedoMap", 0); // 没有贴图，关闭开关 (false)
            }

            // 【新增】：1 号和 2 号插槽：IBL 贴图
            if (m_IrradianceMap && m_PrefilterMap) {
                m_IrradianceMap->Bind(1); // 绑定到 1 号槽
                m_Shader->SetUniform1i("u_IrradianceMap", 1);

                m_PrefilterMap->Bind(2);  // 绑定到 2 号槽
                m_Shader->SetUniform1i("u_PrefilterMap", 2);

                m_Shader->SetUniform1i("u_UseIBL", 1); // 开启 IBL 光照！
            } else {
                m_Shader->SetUniform1i("u_UseIBL", 0);
            }
        }
    };
}
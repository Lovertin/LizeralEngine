#pragma once

#include "Material.h"
#include "runtime/core/math/vector3.h"
#include "runtime/function/res_type/texture/Texture2D.h"
#include "runtime/function/res_type/texture/TextureCube.h"
#include "runtime/core/meta/reflection/reflection.h"
#include "runtime/resource/resourceManager/resourceManager.h" 

namespace Lizeral {
    
    REFLECTION_TYPE(PBRMaterial)
    CLASS(PBRMaterial : public Material, WhiteListFields) {
        REFLECTION_BODY(PBRMaterial)
    public:

        META(Enable)
        Vector3 m_Albedo { 1.0f, 1.0f, 1.0f };
        
        META(Enable)
        float   m_Metallic { 0.0f };          
        
        META(Enable)
        float   m_Roughness { 0.5f };         
        
        META(Enable)
        float   m_AO { 1.0f };                

        META(Enable)
        std::string m_AlbedoMapPath = "";
        
        META(Enable)
        std::string m_IrradianceMapPath = "";
        
        META(Enable)
        std::string m_PrefilterMapPath = "";


        std::shared_ptr<Texture2D> m_AlbedoMap;
        std::shared_ptr<TextureCube> m_IrradianceMap;
        std::shared_ptr<TextureCube> m_PrefilterMap;

        PBRMaterial() = default;
        PBRMaterial(std::shared_ptr<Shader> shader) : Material(shader) {}

        void LoadTextures() {
            auto& rm = ResourceManager::GetInstance();
            if (!m_AlbedoMapPath.empty()) {
                m_AlbedoMap = rm.Load<Texture2D>(m_AlbedoMapPath);
            } else {
                m_AlbedoMap = nullptr; 
            }

            if (!m_IrradianceMapPath.empty()) {
                m_IrradianceMap = rm.Load<TextureCube>(m_IrradianceMapPath);
            } else {
                m_IrradianceMap = nullptr;
            }

            if (!m_PrefilterMapPath.empty()) {
                m_PrefilterMap = rm.Load<TextureCube>(m_PrefilterMapPath);
            } else {
                m_PrefilterMap = nullptr;
            }
        }

        void BindAndApply() override {
            if (!m_Shader) return;
            m_Shader->Bind();

            m_Shader->SetUniformVector3f("u_Albedo", m_Albedo);
            m_Shader->SetUniform1f("u_Metallic", m_Metallic);
            m_Shader->SetUniform1f("u_Roughness", m_Roughness);
            m_Shader->SetUniform1f("u_AO", m_AO);

            if (m_AlbedoMap) {
                m_AlbedoMap->Bind(0); 
                m_Shader->SetUniform1i("u_AlbedoMap", 0); 
                m_Shader->SetUniform1i("u_UseAlbedoMap", 1); 
            } else {
                m_Shader->SetUniform1i("u_UseAlbedoMap", 0); 
            }

            if (m_IrradianceMap && m_PrefilterMap) {
                m_IrradianceMap->Bind(1); 
                m_Shader->SetUniform1i("u_IrradianceMap", 1);
                m_PrefilterMap->Bind(2);  
                m_Shader->SetUniform1i("u_PrefilterMap", 2);
                m_Shader->SetUniform1i("u_UseIBL", 1); 
            } else {
                m_Shader->SetUniform1i("u_UseIBL", 0);
            }
        }
    };
}
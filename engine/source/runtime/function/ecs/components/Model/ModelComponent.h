#pragma once
#include <memory>
#include <string>
#include <vector>
#include "runtime/function/ecs/components/component.h"
#include "runtime/function/res_type/Model/Model.h"
#include "runtime/function/res_type/Material/Material.h"
#include "runtime/resource/resourceManager/resourceManager.h"
#include "runtime/function/res_type/Material/PBRMaterial.h"

namespace Lizeral {
    REFLECTION_TYPE(ModelComponent)
    CLASS(ModelComponent : public Component, WhiteListFields) {
        REFLECTION_BODY(ModelComponent)
    public:

        META(Enable,UI:Address)
        std::string m_ModelPath;

        META(Enable,UI:Address)
        std::vector<std::string> m_OverrideMaterialPaths; 

        META(Enable)
        bool m_UseGlobalMaterial = true;

        std::shared_ptr<Model> m_Model;
        std::string m_LastLoadedModelPath="";
        std::vector<std::shared_ptr<Material>> m_OverrideMaterials; 

        ModelComponent() = default;

        void LoadResources() {
            auto& rm = ResourceManager::GetInstance();

            if (!m_ModelPath.empty()) {
                m_Model = rm.Load<Model>(m_ModelPath);

                if (m_Model) {
                    static std::shared_ptr<Shader> s_defaultPbrShader = nullptr;
                    if (!s_defaultPbrShader) {
                        s_defaultPbrShader = std::make_shared<Shader>(
                            "C:\\Lizeral Engine\\LizeralEngine0.0.1\\engine\\source\\runtime\\function\\sandbox\\shaderTest\\solid.vert",
                            "C:\\Lizeral Engine\\LizeralEngine0.0.1\\engine\\source\\runtime\\function\\sandbox\\shaderTest\\solid.frag"
                        );
                    }

                    for (auto& mesh : m_Model->GetMeshes()) {
                        auto pbrMat = std::dynamic_pointer_cast<PBRMaterial>(mesh.m_Material);
                        if (pbrMat && pbrMat->GetShader() == nullptr) {
                            pbrMat->SetShader(s_defaultPbrShader);
                            pbrMat->LoadTextures(); 
                        }
                    }
                }
            } else {
                m_Model = nullptr;
            }

            m_OverrideMaterials.clear();
            for (const std::string& matPath : m_OverrideMaterialPaths) {
                if (!matPath.empty()) {
                    // m_OverrideMaterials.push_back(rm.Load<PBRMaterial>(matPath));
                } else {
                    m_OverrideMaterials.push_back(nullptr);
                }
            }

            if (m_Model && m_OverrideMaterials.size() < m_Model->GetMeshes().size()) {
                m_OverrideMaterials.resize(m_Model->GetMeshes().size(), nullptr);
            }

            m_LastLoadedModelPath = m_ModelPath;
        }
    };
}
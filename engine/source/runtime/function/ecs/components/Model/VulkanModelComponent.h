#pragma once 
#include "runtime/core/meta/reflection/reflection.h"
#include "runtime/function/ecs/components/component.h"
#include "runtime/resource/material/materialAsset.h"
#include <algorithm>
#include <cstdint>
#include <string>
#include <vector>

namespace Lizeral{

    struct VulkanTextureOverrideSet {
        std::string albedoTexturePath;
        std::string normalTexturePath;
        std::string ormTexturePath;
        std::string emissiveTexturePath;

        bool Empty() const {
            return albedoTexturePath.empty()
                && normalTexturePath.empty()
                && ormTexturePath.empty()
                && emissiveTexturePath.empty();
        }
    };

    struct VulkanMaterialSlotOverride {
        uint32_t meshAssetIndex {0};
        uint32_t materialSlotIndex {0};
        bool enabled {true};
        std::string materialAssetPath;
        Resource::MaterialInstance materialInstance {};
        VulkanTextureOverrideSet textureOverrides {};
    };

    struct VulkanMeshAssetState {
        uint32_t meshAssetIndex {0};
        bool visible {true};
        bool castShadow {true};
    };
        
    REFLECTION_TYPE(VulkanModelComponent)
    CLASS(VulkanModelComponent : public Component, WhiteListFields){
        REFLECTION_BODY(VulkanModelComponent)
    public:
        const std::string& getModelAssetPath() const { return m_ModelPath; }
        std::string getVulkanModelPath() const { return m_ModelPath; }
        bool HasModelAsset() const { return !m_ModelPath.empty(); }
        
        void setModelAssetPath(const std::string& asset_path) {
            if (m_ModelPath != asset_path) {
                m_ModelPath = asset_path;
                MarkResourceDirty();
            }
        }

        void setVulkanModelPath(const std::string& glb_path) {
            setModelAssetPath(glb_path);
        }

        bool IsLoaded() const { return m_IsLoaded; }
        void SetLoaded(bool state) { m_IsLoaded = state; }
        void ResetLoadState() { m_IsLoaded = false; }

        bool IsVisible() const { return m_IsVisible; }
        void SetVisible(bool visible) {
            if (m_IsVisible != visible) {
                m_IsVisible = visible;
                MarkResourceDirty();
            }
        }

        bool CastShadow() const { return m_CastShadow; }
        void SetCastShadow(bool castShadow) {
            if (m_CastShadow != castShadow) {
                m_CastShadow = castShadow;
                MarkResourceDirty();
            }
        }

        bool ReceiveShadow() const { return m_ReceiveShadow; }
        void SetReceiveShadow(bool receiveShadow) {
            if (m_ReceiveShadow != receiveShadow) {
                m_ReceiveShadow = receiveShadow;
                MarkResourceDirty();
            }
        }

        uint32_t GetResourceRevision() const { return m_ResourceRevision; }

        const std::vector<VulkanMaterialSlotOverride>& GetMaterialOverrides() const { return m_MaterialOverrides; }
        std::vector<VulkanMaterialSlotOverride>& GetMaterialOverrides() { return m_MaterialOverrides; }

        const VulkanMaterialSlotOverride* FindMaterialOverride(uint32_t meshAssetIndex, uint32_t materialSlotIndex = 0) const {
            auto it = std::find_if(
                m_MaterialOverrides.begin(),
                m_MaterialOverrides.end(),
                [meshAssetIndex, materialSlotIndex](const VulkanMaterialSlotOverride& overrideData) {
                    return overrideData.meshAssetIndex == meshAssetIndex
                        && overrideData.materialSlotIndex == materialSlotIndex;
                }
            );
            return it != m_MaterialOverrides.end() ? &(*it) : nullptr;
        }

        VulkanMaterialSlotOverride* FindMaterialOverride(uint32_t meshAssetIndex, uint32_t materialSlotIndex = 0) {
            auto it = std::find_if(
                m_MaterialOverrides.begin(),
                m_MaterialOverrides.end(),
                [meshAssetIndex, materialSlotIndex](const VulkanMaterialSlotOverride& overrideData) {
                    return overrideData.meshAssetIndex == meshAssetIndex
                        && overrideData.materialSlotIndex == materialSlotIndex;
                }
            );
            return it != m_MaterialOverrides.end() ? &(*it) : nullptr;
        }

        VulkanMaterialSlotOverride& UpsertMaterialOverride(uint32_t meshAssetIndex, uint32_t materialSlotIndex = 0) {
            if (VulkanMaterialSlotOverride* existing = FindMaterialOverride(meshAssetIndex, materialSlotIndex)) {
                return *existing;
            }

            VulkanMaterialSlotOverride overrideData {};
            overrideData.meshAssetIndex = meshAssetIndex;
            overrideData.materialSlotIndex = materialSlotIndex;
            m_MaterialOverrides.push_back(std::move(overrideData));
            MarkResourceDirty();
            return m_MaterialOverrides.back();
        }

        bool RemoveMaterialOverride(uint32_t meshAssetIndex, uint32_t materialSlotIndex = 0) {
            const auto oldSize = m_MaterialOverrides.size();
            m_MaterialOverrides.erase(
                std::remove_if(
                    m_MaterialOverrides.begin(),
                    m_MaterialOverrides.end(),
                    [meshAssetIndex, materialSlotIndex](const VulkanMaterialSlotOverride& overrideData) {
                        return overrideData.meshAssetIndex == meshAssetIndex
                            && overrideData.materialSlotIndex == materialSlotIndex;
                    }
                ),
                m_MaterialOverrides.end()
            );

            const bool removed = oldSize != m_MaterialOverrides.size();
            if (removed) {
                MarkResourceDirty();
            }
            return removed;
        }

        void ClearMaterialOverrides() {
            if (!m_MaterialOverrides.empty()) {
                m_MaterialOverrides.clear();
                MarkResourceDirty();
            }
        }

        void NotifyMaterialOverridesChanged() {
            MarkResourceDirty();
        }

        const std::vector<VulkanMeshAssetState>& GetMeshAssetStates() const { return m_MeshAssetStates; }
        std::vector<VulkanMeshAssetState>& GetMeshAssetStates() { return m_MeshAssetStates; }

        bool IsMeshVisible(uint32_t meshAssetIndex) const {
            const VulkanMeshAssetState* state = FindMeshState(meshAssetIndex);
            return state ? state->visible : true;
        }

        void SetMeshVisible(uint32_t meshAssetIndex, bool visible) {
            VulkanMeshAssetState& state = UpsertMeshState(meshAssetIndex);
            if (state.visible != visible) {
                state.visible = visible;
                MarkResourceDirty();
            }
        }

        bool DoesMeshCastShadow(uint32_t meshAssetIndex) const {
            const VulkanMeshAssetState* state = FindMeshState(meshAssetIndex);
            return state ? state->castShadow : true;
        }

        void SetMeshCastShadow(uint32_t meshAssetIndex, bool castShadow) {
            VulkanMeshAssetState& state = UpsertMeshState(meshAssetIndex);
            if (state.castShadow != castShadow) {
                state.castShadow = castShadow;
                MarkResourceDirty();
            }
        }

        void NotifyMeshStatesChanged() {
            MarkResourceDirty();
        }

        void onReflectionUpdated(const char* field_name) {
            if (field_name == nullptr) {
                return;
            }

            const std::string fieldName(field_name);
            if (fieldName == "m_ModelPath" || fieldName == "ModelPath"
                || fieldName == "m_IsVisible" || fieldName == "Visible"
                || fieldName == "m_CastShadow" || fieldName == "CastShadow"
                || fieldName == "m_ReceiveShadow" || fieldName == "ReceiveShadow") {
                MarkResourceDirty();
            }
        }

    private:
        const VulkanMeshAssetState* FindMeshState(uint32_t meshAssetIndex) const {
            auto it = std::find_if(
                m_MeshAssetStates.begin(),
                m_MeshAssetStates.end(),
                [meshAssetIndex](const VulkanMeshAssetState& state) {
                    return state.meshAssetIndex == meshAssetIndex;
                }
            );
            return it != m_MeshAssetStates.end() ? &(*it) : nullptr;
        }

        VulkanMeshAssetState* FindMeshState(uint32_t meshAssetIndex) {
            auto it = std::find_if(
                m_MeshAssetStates.begin(),
                m_MeshAssetStates.end(),
                [meshAssetIndex](const VulkanMeshAssetState& state) {
                    return state.meshAssetIndex == meshAssetIndex;
                }
            );
            return it != m_MeshAssetStates.end() ? &(*it) : nullptr;
        }

        VulkanMeshAssetState& UpsertMeshState(uint32_t meshAssetIndex) {
            if (VulkanMeshAssetState* state = FindMeshState(meshAssetIndex)) {
                return *state;
            }

            VulkanMeshAssetState state {};
            state.meshAssetIndex = meshAssetIndex;
            m_MeshAssetStates.push_back(state);
            return m_MeshAssetStates.back();
        }

        void MarkResourceDirty() {
            ++m_ResourceRevision;
            m_IsLoaded = false;
        }

        META(Enable, UI:Address)
        std::string m_ModelPath;

        META(Enable)
        bool m_IsVisible = true;
        META(Enable)
        bool m_CastShadow = true;
        META(Enable)
        bool m_ReceiveShadow = true;
        bool m_IsLoaded = false; 
        uint32_t m_ResourceRevision = 0;
        std::vector<VulkanMaterialSlotOverride> m_MaterialOverrides;
        std::vector<VulkanMeshAssetState> m_MeshAssetStates;
    };

}

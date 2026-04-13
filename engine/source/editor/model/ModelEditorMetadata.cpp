#include "ModelEditorMetadata.h"

#include "runtime/resource/resourceLoadingUnit.h"

#include <iostream>
#include <unordered_map>

namespace Lizeral {

    namespace {
        std::unordered_map<std::string, EditorModelMetadata> s_metadataCache;
    }

    const EditorModelMetadata* GetModelEditorMetadata(const std::string& modelPath) {
        if (modelPath.empty()) {
            return nullptr;
        }

        auto cached = s_metadataCache.find(modelPath);
        if (cached != s_metadataCache.end()) {
            return &cached->second;
        }

        EditorModelMetadata metadata;
        Resource::ResourceLoadingUnit loader;
        Resource::RuntimeModelData runtimeModel;
        std::string error;

        if (loader.LoadModel(modelPath, &runtimeModel, 0, &error)) {
            metadata.valid = true;
            metadata.modelName = runtimeModel.modelName;
            metadata.materialNames.reserve(runtimeModel.materialAssets.size());

            for (size_t materialIndex = 0; materialIndex < runtimeModel.materialAssets.size(); ++materialIndex) {
                std::string materialName = runtimeModel.materialAssets[materialIndex].name;
                if (materialName.empty()) {
                    materialName = "Material " + std::to_string(materialIndex);
                }
                metadata.materialNames.push_back(materialName);
            }

            metadata.meshEntries.reserve(runtimeModel.meshAssets.size());
            for (size_t meshIndex = 0; meshIndex < runtimeModel.meshAssets.size(); ++meshIndex) {
                const auto& meshAsset = runtimeModel.meshAssets[meshIndex];

                EditorModelMeshEntry entry;
                entry.meshAssetIndex = static_cast<uint32_t>(meshIndex);
                entry.nodeIndex = meshAsset.nodeIndex;
                entry.materialIndex = meshAsset.materialIndex;
                entry.meshName = meshAsset.name.empty() ? ("SubMesh " + std::to_string(meshIndex)) : meshAsset.name;
                if (meshAsset.nodeIndex < runtimeModel.nodes.size()) {
                    entry.nodeName = runtimeModel.nodes[meshAsset.nodeIndex].name;
                }
                entry.materialName = entry.materialIndex < metadata.materialNames.size()
                    ? metadata.materialNames[entry.materialIndex]
                    : ("Material " + std::to_string(entry.materialIndex));

                entry.displayName = entry.meshName;
                if (!entry.nodeName.empty()) {
                    entry.displayName += " [" + entry.nodeName + "]";
                }
                entry.displayName += " -> " + entry.materialName;
                metadata.meshEntries.push_back(std::move(entry));
            }
        } else {
            std::cerr << "[ModelEditorMetadata] WARNING: Failed to inspect model '" << modelPath
                      << "'. Reason: " << error << std::endl;
        }

        auto [it, inserted] = s_metadataCache.emplace(modelPath, std::move(metadata));
        (void)inserted;
        return &it->second;
    }

    void ClearModelEditorMetadataCache(const std::string& modelPath) {
        if (modelPath.empty()) {
            s_metadataCache.clear();
        } else {
            s_metadataCache.erase(modelPath);
        }
    }

} // namespace Lizeral

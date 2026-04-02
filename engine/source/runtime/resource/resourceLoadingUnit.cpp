#include "runtime/resource/resourceLoadingUnit.h"

#include <iostream>

namespace Lizeral::Resource {

    bool ResourceLoadingUnit::LoadModel(const std::string& sourcePath, RuntimeModelData* outModel, uint32_t globalTextureOffset, std::string* outError) const {
        ImportedModelData importedModel;
        std::string importError;
        if (!m_importer.ImportModel(sourcePath, &importedModel, &importError)) {
            if (outError) {
                *outError = importError.empty() ? "Importer failed with unknown error." : importError;
            }
            return false;
        }

        MeshletCookOptions cookOptions{};
        cookOptions.globalTextureOffset = globalTextureOffset;
        cookOptions.maxVerticesPerMeshlet = 64;
        cookOptions.maxTrianglesPerMeshlet = 124;
        cookOptions.ensureFallbackTextureForTexturelessModel = true;

        std::string cookError;
        if (!m_cooker.Cook(importedModel, cookOptions, outModel, &cookError)) {
            if (outError) {
                *outError = cookError.empty() ? "Cooker failed with unknown error." : cookError;
            }
            return false;
        }

        std::cout << "[ResourceLoadingUnit] Model loaded successfully: " << sourcePath << std::endl;
        return true;
    }

} // namespace Lizeral::Resource


#pragma once

#include "runtime/resource/modelResourceData.h"
#include "runtime/resource/cooker/meshletModelCooker.h"
#include "runtime/resource/importer/assimpModelImporter.h"

#include <string>

namespace Lizeral::Resource {

    class ResourceLoadingUnit {
    public:
        bool LoadModel(const std::string& sourcePath, RuntimeModelData* outModel, uint32_t globalTextureOffset = 0, std::string* outError = nullptr) const;

    private:
        AssimpModelImporter m_importer;
        MeshletModelCooker m_cooker;
    };

} // namespace Lizeral::Resource


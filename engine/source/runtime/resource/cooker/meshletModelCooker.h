#pragma once

#include "runtime/resource/modelResourceData.h"

#include <string>

namespace Lizeral::Resource {

    class MeshletModelCooker {
    public:
        bool Cook(const ImportedModelData& importedModel, const MeshletCookOptions& options, RuntimeModelData* outModel, std::string* outError = nullptr) const;
    };

} // namespace Lizeral::Resource


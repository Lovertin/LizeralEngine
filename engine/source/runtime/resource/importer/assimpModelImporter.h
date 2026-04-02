#pragma once

#include "runtime/resource/modelResourceData.h"

#include <string>

namespace Lizeral::Resource {

    class AssimpModelImporter {
    public:
        bool ImportModel(const std::string& sourcePath, ImportedModelData* outModel, std::string* outError = nullptr) const;
    };

} // namespace Lizeral::Resource


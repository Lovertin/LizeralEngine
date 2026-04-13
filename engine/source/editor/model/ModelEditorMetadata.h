#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace Lizeral {

    struct EditorModelMeshEntry {
        uint32_t meshAssetIndex = 0;
        uint32_t nodeIndex = 0;
        uint32_t materialIndex = 0;
        std::string meshName;
        std::string nodeName;
        std::string materialName;
        std::string displayName;
    };

    struct EditorModelMetadata {
        bool valid = false;
        std::string modelName;
        std::vector<EditorModelMeshEntry> meshEntries;
        std::vector<std::string> materialNames;
    };

    const EditorModelMetadata* GetModelEditorMetadata(const std::string& modelPath);
    void ClearModelEditorMetadataCache(const std::string& modelPath = {});

} // namespace Lizeral

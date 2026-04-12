#pragma once
#include <string>
#include "runtime/core/meta/reflection/reflection.h"

namespace Lizeral {
    class ResourceAsset {
    public:
        virtual ~ResourceAsset() = default;

        const std::string& GetPath() const { return m_path; }
        void SetPath(const std::string& path) { m_path = path; }

        virtual bool LoadFromFile(const std::string& path) = 0;

    protected:
        std::string m_path {""};
    };
}

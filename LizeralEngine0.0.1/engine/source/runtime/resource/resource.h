#pragma once
#include <string>
#include "runtime/core/meta/reflection/reflection.h"

namespace Lizeral {

    REFLECTION_TYPE(Resource)
    CLASS(Resource, WhiteListFields) {
        REFLECTION_BODY(Resource)

    public:
        virtual ~Resource() = default;

        const std::string& GetPath() const { return m_path; }
        void SetPath(const std::string& path) { m_path = path; }

        // 【核心契约】：所有继承我的资源，都必须实现这个函数！
        // ResourceManager 会在不知道你是谁的情况下，无脑调用这个函数。
        virtual bool LoadFromFile(const std::string& path) = 0;

    protected:
        META(Enable)
        std::string m_path {""};
    };
}
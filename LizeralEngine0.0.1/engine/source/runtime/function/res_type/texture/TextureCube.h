#pragma once

#include <string>
#include <vector>
#include <glad/glad.h>
#include "runtime/resource/resource.h" 
#include "runtime/core/meta/reflection/reflection.h"


namespace Lizeral {
    REFLECTION_TYPE(TextureCube)
    CLASS(TextureCube : public Resource,WhiteListFields) {
        REFLECTION_BODY(TextureCube)
    private:
        GLuint m_RendererID;
        int m_Width;
        int m_Height;
        int m_Channels;

    public:
        TextureCube();
        virtual ~TextureCube();

        // 这里传入的 path 是基础路径，例如 "asset/skybox/skybox_irradiance"
        // 内部会自动拼接 "_X+.hdr" 等后缀来加载 6 张图片
        bool LoadFromFile(const std::string& basePath) override;

        void Bind(uint32_t slot = 0) const;
        void Unbind() const;

        GLuint GetRendererID() const { return m_RendererID; }
    };

}
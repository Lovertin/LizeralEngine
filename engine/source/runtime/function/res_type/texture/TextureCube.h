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
        ~TextureCube();

        bool LoadFromFile(const std::string& basePath) override;

        void Bind(uint32_t slot = 0) const;
        void Unbind() const;

        GLuint GetRendererID() const { return m_RendererID; }
    };

}
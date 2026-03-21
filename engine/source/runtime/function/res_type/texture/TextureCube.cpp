#include "TextureCube.h"
#include <iostream>

#include "runtime/core/vendor/stb_image/stb_image.h"

namespace Lizeral {

    TextureCube::TextureCube() : m_RendererID(0), m_Width(0), m_Height(0), m_Channels(0) {
    }

    TextureCube::~TextureCube() {
        if (m_RendererID != 0) {
            glDeleteTextures(1, &m_RendererID);
        }
    }

    bool TextureCube::LoadFromFile(const std::string& basePath) {
        glGenTextures(1, &m_RendererID);
        // GL_TEXTURE_CUBE_MAP
        glBindTexture(GL_TEXTURE_CUBE_MAP, m_RendererID);

        std::vector<std::string> suffixes = {
            "_X+.hdr", "_X-.hdr",
            "_Z+.hdr", "_Z-.hdr",
            "_Y+.hdr", "_Y-.hdr"
        };

        stbi_set_flip_vertically_on_load(false);

        for (unsigned int i = 0; i < 6; i++) {
            std::string fullPath = basePath + suffixes[i];
            
            float* data = stbi_loadf(fullPath.c_str(), &m_Width, &m_Height, &m_Channels, 0);
            
            if (data) {
                GLenum internalFormat = (m_Channels == 4) ? GL_RGBA16F : GL_RGB16F;
                GLenum dataFormat = (m_Channels == 4) ? GL_RGBA : GL_RGB;

                glTexImage2D(
                    GL_TEXTURE_CUBE_MAP_POSITIVE_X + i, 
                    0, internalFormat, m_Width, m_Height, 0, dataFormat, GL_FLOAT, data 
                );
                
                stbi_image_free(data);
            } else {
                std::cerr << "[TextureCube] ERROR: Failed to load cubemap face: " << fullPath << std::endl;
                return false;
            }
        }
        glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        
        glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);

        std::cout << "[TextureCube] Successfully loaded HDRI cubemap from base path: " << basePath << std::endl;
        return true;
    }

    void TextureCube::Bind(uint32_t slot) const {
        glActiveTexture(GL_TEXTURE0 + slot);
        glBindTexture(GL_TEXTURE_CUBE_MAP, m_RendererID);
    }

    void TextureCube::Unbind() const {
        glBindTexture(GL_TEXTURE_CUBE_MAP, 0);
    }

}
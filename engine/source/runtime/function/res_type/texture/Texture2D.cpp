#include "Texture2D.h"
#include <iostream>

#include "runtime/core/vendor/stb_image/stb_image.h" 

namespace Lizeral {

    Texture2D::Texture2D()
        : m_RendererID(0), m_Width(0), m_Height(0), m_Channels(0) {
    }

    Texture2D::~Texture2D() {
        if (m_RendererID != 0) {
            glDeleteTextures(1, &m_RendererID);
        }
    }

    bool Texture2D::LoadFromFile(const std::string& path) {
        m_Path = path;

        stbi_set_flip_vertically_on_load(true);

        glGenTextures(1, &m_RendererID);
        glBindTexture(GL_TEXTURE_2D, m_RendererID);

        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

        bool isHDR = stbi_is_hdr(path.c_str());

        if (isHDR) {
            float* data = stbi_loadf(path.c_str(), &m_Width, &m_Height, &m_Channels, 0);
            if (!data) {
                std::cerr << "[Texture2D] Failed to load HDR texture: " << path << std::endl;
                return false;
            }

            GLenum internalFormat = (m_Channels == 4) ? GL_RGBA16F : GL_RGB16F;
            GLenum dataFormat = (m_Channels == 4) ? GL_RGBA : GL_RGB;

            glTexImage2D(GL_TEXTURE_2D, 0, internalFormat, m_Width, m_Height, 0, dataFormat, GL_FLOAT, data);
            glGenerateMipmap(GL_TEXTURE_2D);
            stbi_image_free(data);

        } else {
            stbi_uc* data = stbi_load(path.c_str(), &m_Width, &m_Height, &m_Channels, 0);
            if (!data) {
                std::cerr << "[Texture2D] Failed to load texture: " << path << std::endl;
                return false;
            }

            GLenum internalFormat = 0, dataFormat = 0;
            if (m_Channels == 4) {
                internalFormat = GL_RGBA8; dataFormat = GL_RGBA;
            } else if (m_Channels == 3) {
                internalFormat = GL_RGB8; dataFormat = GL_RGB;
            } else if (m_Channels == 1) {
                internalFormat = GL_R8; dataFormat = GL_RED;
            }

            glTexImage2D(GL_TEXTURE_2D, 0, internalFormat, m_Width, m_Height, 0, dataFormat, GL_UNSIGNED_BYTE, data);
            glGenerateMipmap(GL_TEXTURE_2D);
            stbi_image_free(data);
        }

        std::cout << "[Texture] Loaded successfully: " << path << " (" << m_Width << "x" << m_Height << (isHDR ? " HDR" : " LDR") << ")" << std::endl;
        return true;
    }

    bool Texture2D::LoadFromMemory(const unsigned char* buffer, int length) {
        glGenTextures(1, &m_RendererID);
        glBindTexture(GL_TEXTURE_2D, m_RendererID);

        stbi_set_flip_vertically_on_load(false); 

        int width, height, channels;
        unsigned char* data = stbi_load_from_memory(buffer, length, &width, &height, &channels, 0);

        if (data) {
            GLenum format = (channels == 4) ? GL_RGBA : GL_RGB;
            glTexImage2D(GL_TEXTURE_2D, 0, format, width, height, 0, format, GL_UNSIGNED_BYTE, data);
            glGenerateMipmap(GL_TEXTURE_2D);

            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

            stbi_image_free(data);
            return true;
        } else {
            std::cerr << "[Texture2D] ERROR: Failed to load texture from memory!" << std::endl;
            return false;
        }
    }

    void Texture2D::Bind(uint32_t slot) const {
        glActiveTexture(GL_TEXTURE0 + slot);
        glBindTexture(GL_TEXTURE_2D, m_RendererID);
    }

    void Texture2D::Unbind() const {
        glBindTexture(GL_TEXTURE_2D, 0);
    }

}
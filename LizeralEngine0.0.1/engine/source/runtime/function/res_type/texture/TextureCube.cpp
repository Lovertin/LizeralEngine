#include "TextureCube.h"
#include <iostream>

// 引入 stb_image
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
        // 注意：这里绑定的目标变成了 GL_TEXTURE_CUBE_MAP
        glBindTexture(GL_TEXTURE_CUBE_MAP, m_RendererID);

        // OpenGL 的天空盒面固定顺序：右(X+)、左(X-)、上(Y+)、下(Y-)、后(Z+)、前(Z-)
        // 这个顺序必须严格和 GL_TEXTURE_CUBE_MAP_POSITIVE_X 等枚举的递增顺序对应
        std::vector<std::string> suffixes = {
            "_X+.hdr", "_X-.hdr",
            "_Z+.hdr", "_Z-.hdr",
            "_Y+.hdr", "_Y-.hdr"
        };

        // HDRI Cubemap 采样坐标系和普通 2D 贴图不同，一般不需要（也不能）在加载时翻转 Y 轴
        stbi_set_flip_vertically_on_load(false);

        for (unsigned int i = 0; i < 6; i++) {
            std::string fullPath = basePath + suffixes[i];
            
            // 【硬核节点】：使用 stbi_loadf 读取浮点数 HDR 数据！
            float* data = stbi_loadf(fullPath.c_str(), &m_Width, &m_Height, &m_Channels, 0);
            
            if (data) {
                // 判断格式，HDR 通常使用 16位 浮点数存入显存 (GL_RGB16F)
                GLenum internalFormat = (m_Channels == 4) ? GL_RGBA16F : GL_RGB16F;
                GLenum dataFormat = (m_Channels == 4) ? GL_RGBA : GL_RGB;

                // 连续绑定 6 个面
                glTexImage2D(
                    GL_TEXTURE_CUBE_MAP_POSITIVE_X + i, 
                    0, internalFormat, m_Width, m_Height, 0, dataFormat, GL_FLOAT, data // 数据类型变成 GL_FLOAT
                );
                
                stbi_image_free(data);
            } else {
                std::cerr << "[TextureCube] ERROR: Failed to load cubemap face: " << fullPath << std::endl;
                return false;
            }
        }

        // 设置采样参数：放大和缩小均使用线性插值
        glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        
        // 【关键】：三轴必须设置为 CLAMP_TO_EDGE，消除天空盒折角处的接缝黑线！
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
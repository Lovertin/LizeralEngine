#include "Texture2D.h"
#include <iostream>

// 引入刚刚配置好的 stb_image
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

        // 【大坑预警】：OpenGL 的纹理坐标 (0,0) 在左下角，而大多数图片格式的 (0,0) 在左上角！
        // 必须告诉 stbi 在加载时把图片上下翻转，否则贴图会是倒着的！
        stbi_set_flip_vertically_on_load(true);

        // 1. 从硬盘读取像素数据到内存
        stbi_uc* data = stbi_load(path.c_str(), &m_Width, &m_Height, &m_Channels, 0);
        
        if (!data) {
            std::cerr << "Failed to load texture: " << path << std::endl;
            return false;
        }

        // 判断图片是 RGB 还是 RGBA (有没有透明通道)
        GLenum internalFormat = 0, dataFormat = 0;
        if (m_Channels == 4) {
            internalFormat = GL_RGBA8;
            dataFormat = GL_RGBA;
        } else if (m_Channels == 3) {
            internalFormat = GL_RGB8;
            dataFormat = GL_RGB;
        } else if (m_Channels == 1) { // 比如单独的粗糙度灰度图
            internalFormat = GL_R8;
            dataFormat = GL_RED;
        } else {
            std::cerr << "Unsupported texture channels: " << m_Channels << " in " << path << std::endl;
            stbi_image_free(data);
            return false;
        }

        // 2. 将数据送入显存 (VRAM)
        glGenTextures(1, &m_RendererID);
        glBindTexture(GL_TEXTURE_2D, m_RendererID);

        // 设置纹理环绕方式 (超出 UV 坐标 0~1 的部分怎么处理) - PBR 通常用重复
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        
        // 设置纹理过滤方式 (图片被放大或缩小时怎么采样)
        // 缩小使用 Mipmap 线性插值，放大使用线性插值
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

        // 上传像素数据给 GPU
        glTexImage2D(GL_TEXTURE_2D, 0, internalFormat, m_Width, m_Height, 0, dataFormat, GL_UNSIGNED_BYTE, data);
        
        // 自动生成 Mipmap（多级渐远纹理，极大提升远处物体的渲染性能和抗锯齿）
        glGenerateMipmap(GL_TEXTURE_2D);

        // 3. 清理 CPU 内存 (数据已经进显卡了，RAM 里的可以扔掉了)
        stbi_image_free(data);

        std::cout << "[Texture] Loaded successfully: " << path << " (" << m_Width << "x" << m_Height << ")" << std::endl;
        return true;
    }

    void Texture2D::Bind(uint32_t slot) const {
        // 激活特定的纹理插槽 (OpenGL 至少支持 16 个插槽，GL_TEXTURE0, GL_TEXTURE1...)
        glActiveTexture(GL_TEXTURE0 + slot);
        glBindTexture(GL_TEXTURE_2D, m_RendererID);
    }

    void Texture2D::Unbind() const {
        glBindTexture(GL_TEXTURE_2D, 0);
    }

}
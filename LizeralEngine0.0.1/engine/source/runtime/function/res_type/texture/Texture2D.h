#pragma once

#include <string>
#include <glad/glad.h>
// 如果你的 ResourceManager 强制要求基类，这里可以继承 public Resource
#include "runtime/core/meta/reflection/reflection.h"

#include "runtime/resource/resource.h"

namespace Lizeral {

    REFLECTION_TYPE(Texture2D)
    CLASS(Texture2D:public Resource,WhiteListFields) {
        REFLECTION_BODY(Texture2D)
    private:
        GLuint m_RendererID;  // OpenGL 分配的纹理 ID
        std::string m_Path;   // 文件路径，方便调试
        int m_Width;
        int m_Height;
        int m_Channels;       // 颜色通道数 (RGB=3, RGBA=4)

    public:
        Texture2D();
        virtual ~Texture2D();

        // 核心：供给 ResourceManager 调用的加载接口
        virtual bool LoadFromFile(const std::string& path) override;
        
        bool LoadFromMemory(const unsigned char* buffer, int length);

        // 绑定该贴图到显卡的指定纹理槽位 (Slot)，PBR 通常需要同时绑定多张贴图
        void Bind(uint32_t slot = 0) const;
        void Unbind() const;

        // Getters
        int GetWidth() const { return m_Width; }
        int GetHeight() const { return m_Height; }
        GLuint GetRendererID() const { return m_RendererID; }
    };

}
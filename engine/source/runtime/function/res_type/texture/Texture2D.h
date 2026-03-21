#pragma once

#include <string>
#include <glad/glad.h>
#include "runtime/core/meta/reflection/reflection.h"

#include "runtime/resource/resource.h"

namespace Lizeral {

    REFLECTION_TYPE(Texture2D)
    CLASS(Texture2D:public Resource,WhiteListFields) {
        REFLECTION_BODY(Texture2D)
    private:
        GLuint m_RendererID;  

        META(Enable)
        std::string m_Path;  

        int m_Width;
        int m_Height;
        int m_Channels;      

    public:
        Texture2D();
        ~Texture2D();

        virtual bool LoadFromFile(const std::string& path) override;
        
        bool LoadFromMemory(const unsigned char* buffer, int length);

        void Bind(uint32_t slot = 0) const;
        void Unbind() const;

        // Getters
        int GetWidth() const { return m_Width; }
        int GetHeight() const { return m_Height; }
        GLuint GetRendererID() const { return m_RendererID; }
    };

}
#pragma once
#include <string>
#include "runtime/core/math/vector2.h"
#include "runtime/core/math/vector3.h"
#include "runtime/core/math/vector4.h"
#include "runtime/core/math/matrix3.h"
#include "runtime/core/math/matrix4.h"
// #include "runtime/core/math/math.h"
#include <unordered_map>
#include <glad/glad.h>

namespace Lizeral{
    class Shader{

    public:
        Shader(const std::string& vertexPath, const std::string& fragmentPath);
        virtual ~Shader();

        Shader(const Shader&) = delete;
        Shader& operator=(const Shader&) = delete;
        Shader(Shader&& other) noexcept;
        Shader& operator=(Shader&& other) noexcept;

        // 绑定/解绑着色器
        void Bind() const;
        void Unbind() const;

        void SetUniform1i(const std::string& name, int value);
        void SetUniform1f(const std::string& name, float value);
        void SetUniform2f(const std::string& name, float v0, float v1);
        void SetUniform3f(const std::string& name, float v0, float v1, float v2);
        void SetUniform4f(const std::string& name, float v0, float v1, float v2, float v3);

        void SetUniformMat3f(const std::string& name,const Matrix3x3& mat3);
        void SetUniformMat4f(const std::string& name,const Matrix4x4& mat4);
        void SetUniformVector2f(const std::string& name,const Vector2& vec2);
        void SetUniformVector3f(const std::string& name,const Vector3& vec3);
        void SetUniformVector4f(const std::string& name,const Vector4& vec4);

    private:

        std::string ReadFile(const std::string& filepath);
        unsigned int CompileShader(unsigned int type, const std::string& source);
        unsigned int CreateProgram(const std::string& vertexSource, const std::string& fragmentSource);
        int GetUniformLocation(const std::string& name);

        void CheckShaderErrors(unsigned int shader, unsigned int type);
        void CheckProgramErrors(unsigned int program);

    private:
        unsigned int m_rendererID;  // 着色器程序ID
        std::unordered_map<std::string, int> m_UniformLocationCache;  // uniform变量位置缓存
    };
}
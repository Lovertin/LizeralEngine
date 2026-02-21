#include "shader.h"
#include<iostream>
#include <fstream>
#include <sstream>
namespace Lizeral{
    Shader::Shader(const std::string& vertexPath, const std::string& fragmentPath) {
        // 读取着色器源码
        std::string vertexSource = ReadFile(vertexPath);
        std::string fragmentSource = ReadFile(fragmentPath);
        
        // 创建并链接着色器程序
        m_rendererID = CreateProgram(vertexSource, fragmentSource);
    }

    Shader::~Shader() {
        glDeleteProgram(m_rendererID);
    }   

    Shader::Shader(Shader&& other) noexcept 
    : m_rendererID(other.m_rendererID), 
      m_UniformLocationCache(std::move(other.m_UniformLocationCache)) {
        other.m_rendererID = 0;
    }

    Shader& Shader::operator=(Shader&& other) noexcept {
        if (this != &other) {
            glDeleteProgram(m_rendererID);
            m_rendererID = other.m_rendererID;
            m_UniformLocationCache = std::move(other.m_UniformLocationCache);
            other.m_rendererID = 0;
        }
        return *this;
    }

    void Shader::Bind() const {
        glUseProgram(m_rendererID);
    }

    void Shader::Unbind() const {
        glUseProgram(0);
    }

    // Uniform设置实现
    void Shader::SetUniform1i(const std::string& name, int value) {
        glUniform1i(GetUniformLocation(name), value);
    }

    void Shader::SetUniform1f(const std::string& name, float value) {
        glUniform1f(GetUniformLocation(name), value);
    }

    void Shader::SetUniform2f(const std::string& name, float v0, float v1) {
        glUniform2f(GetUniformLocation(name), v0, v1);
    }

    void Shader::SetUniform3f(const std::string& name, float v0, float v1, float v2) {
        glUniform3f(GetUniformLocation(name), v0, v1, v2);
    }

    void Shader::SetUniform4f(const std::string& name, float v0, float v1, float v2, float v3) {
        glUniform4f(GetUniformLocation(name), v0, v1, v2, v3);
    }

    void Shader::SetUniformMat3f(const std::string& name,const Matrix3x3& mat3){
        glUniformMatrix3fv(GetUniformLocation(name),1,GL_TRUE,&mat3[0][0]);
    }

    void Shader::SetUniformMat4f(const std::string& name,const Matrix4x4& mat4){
        glUniformMatrix4fv(GetUniformLocation(name),1,GL_TRUE,&mat4[0][0]);
    }

    void Shader::SetUniformVector2f(const std::string& name,const Vector2& vec2){
        glUniform2f(GetUniformLocation(name), vec2.x, vec2.y);
    }

    void Shader::SetUniformVector3f(const std::string& name,const Vector3& vec3){
        glUniform3f(GetUniformLocation(name),vec3.x,vec3.y,vec3.z);
    }

    void Shader::SetUniformVector4f(const std::string& name,const Vector4& vec4){
        glUniform4f(GetUniformLocation(name),vec4.x,vec4.y,vec4.z,vec4.w);
    }

    //辅助函数实现
    std::string Shader::ReadFile(const std::string& filepath) {
        std::ifstream file(filepath);
        std::stringstream stream;
        
        if (!file.is_open()) {
            std::cerr << "Failed to open file: " << filepath << std::endl;
            return "";
        }
        
        stream << file.rdbuf();
        return stream.str();
    }

    unsigned int Shader::CompileShader(unsigned int type, const std::string& source) {
        unsigned int shader = glCreateShader(type);
        const char* src = source.c_str();
        glShaderSource(shader, 1, &src, nullptr);
        glCompileShader(shader);
        
        // 检查编译错误
        CheckShaderErrors(shader, type);
        
        return shader;
    }

    unsigned int Shader::CreateProgram(const std::string& vertexSource, const std::string& fragmentSource) {
        // 编译着色器
        unsigned int vertexShader = CompileShader(GL_VERTEX_SHADER, vertexSource);
        unsigned int fragmentShader = CompileShader(GL_FRAGMENT_SHADER, fragmentSource);
        
        // 创建程序
        unsigned int program = glCreateProgram();
        glAttachShader(program, vertexShader);
        glAttachShader(program, fragmentShader);
        glLinkProgram(program);
        
        // 检查链接错误
        CheckProgramErrors(program);
        
        // 清理着色器对象
        glDeleteShader(vertexShader);
        glDeleteShader(fragmentShader);
        
        return program;
    }

    int Shader::GetUniformLocation(const std::string& name) {
        // 查找缓存
        auto it = m_UniformLocationCache.find(name);
        if (it != m_UniformLocationCache.end()) {
            return it->second;
        }
        
        // 如果不在缓存中，获取位置并存储
        int location = glGetUniformLocation(m_rendererID, name.c_str());
        if (location == -1) {
            std::cerr << "Warning: uniform '" << name << "' doesn't exist!" << std::endl;
        }
        
        m_UniformLocationCache[name] = location;
        return location;
    }

    void Shader::CheckShaderErrors(unsigned int shader, unsigned int type) {
        int success;
        char infoLog[512];
        
        glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
        if (!success) {
            glGetShaderInfoLog(shader, 512, NULL, infoLog);
            std::string typeName = (type == GL_VERTEX_SHADER) ? "VERTEX" : "FRAGMENT";
            std::cerr << "ERROR::SHADER::" << typeName << "::COMPILATION_FAILED\n" << infoLog << std::endl;
        }
    }

    void Shader::CheckProgramErrors(unsigned int program) {
        int success;
        char infoLog[512];
        
        glGetProgramiv(program, GL_LINK_STATUS, &success);
        if (!success) {
            glGetProgramInfoLog(program, 512, NULL, infoLog);
            std::cerr << "ERROR::PROGRAM::LINKING_FAILED\n" << infoLog << std::endl;
        }
    }
}


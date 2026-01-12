#pragma once

#include "runtime/core/meta/reflection/reflection.h"
#include "runtime/core/math/math.h"
#include <glm/glm.hpp>
#include <cassert>

namespace Lizeral {
    REFLECTION_TYPE(Vector3)
    CLASS(Vector3, Fields) {
        REFLECTION_BODY(Vector3)
    
    private:
        glm::vec3 m_data;
    
    public:
        // 构造函数
        Vector3() : m_data(0.0f, 0.0f, 0.0f) {}
        Vector3(float x, float y, float z) : m_data(x, y, z) {}
        explicit Vector3(const glm::vec3& v) : m_data(v) {}
        explicit Vector3(float scalar) : m_data(scalar, scalar, scalar) {}
        
        // 字段访问（用于反射）
        float getX() const { return m_data.x; }
        void setX(float value) { m_data.x = value; }
        
        float getY() const { return m_data.y; }
        void setY(float value) { m_data.y = value; }
        
        float getZ() const { return m_data.z; }
        void setZ(float value) { m_data.z = value; }
        
        // 基本运算符
        Vector3 operator+(const Vector3& rhs) const {
            return Vector3(m_data + rhs.m_data);
        }
        
        Vector3 operator-(const Vector3& rhs) const {
            return Vector3(m_data - rhs.m_data);
        }
        
        Vector3 operator*(float scalar) const {
            return Vector3(m_data * scalar);
        }
        
        Vector3 operator*(const Vector3& rhs) const {
            return Vector3(m_data * rhs.m_data);
        }
        
        Vector3 operator/(float scalar) const {
            assert(scalar != 0.0f);
            return Vector3(m_data / scalar);
        }
        
        Vector3 operator/(const Vector3& rhs) const {
            return Vector3(m_data / rhs.m_data);
        }
        
        const Vector3& operator+() const { return *this; }
        Vector3 operator-() const { return Vector3(-m_data); }
        
        // 复合赋值
        Vector3& operator+=(const Vector3& rhs) {
            m_data += rhs.m_data;
            return *this;
        }
        
        Vector3& operator-=(const Vector3& rhs) {
            m_data -= rhs.m_data;
            return *this;
        }
        
        Vector3& operator*=(float scalar) {
            m_data *= scalar;
            return *this;
        }
        
        Vector3& operator/=(float scalar) {
            assert(scalar != 0.0f);
            m_data /= scalar;
            return *this;
        }
        
        // 友元函数
        friend Vector3 operator*(float scalar, const Vector3& rhs) {
            return rhs * scalar;
        }
        
        // 比较运算符
        bool operator==(const Vector3& rhs) const {
            return m_data == rhs.m_data;
        }
        
        bool operator!=(const Vector3& rhs) const {
            return m_data != rhs.m_data;
        }
        
        // 下标访问
        float operator[](size_t i) const {
            assert(i < 3);
            return (i == 0 ? m_data.x : (i == 1 ? m_data.y : m_data.z));
        }
        
        float& operator[](size_t i) {
            assert(i < 3);
            return (i == 0 ? m_data.x : (i == 1 ? m_data.y : m_data.z));
        }
        
        // 向量运算
        float length() const {
            return glm::length(m_data);
        }
        
        float squaredLength() const {
            return glm::dot(m_data, m_data);
        }
        
        float distance(const Vector3& rhs) const {
            return glm::distance(m_data, rhs.m_data);
        }
        
        float squaredDistance(const Vector3& rhs) const {
            glm::vec3 diff = m_data - rhs.m_data;
            return glm::dot(diff, diff);
        }
        
        float dot(const Vector3& vec) const {
            return glm::dot(m_data, vec.m_data);
        }
        
        Vector3 cross(const Vector3& vec) const {
            return Vector3(glm::cross(m_data, vec.m_data));
        }
        
        // 归一化
        float normalise() {
            float len = glm::length(m_data);
            if (len > 0.0f) {
                float inv_length = 1.0f / len;
                m_data *= inv_length;
            }
            return len;
        }
        
        Vector3 normalisedCopy() const {
            Vector3 ret = *this;
            ret.normalise();
            return ret;
        }
        
        // 检查是否为零向量
        bool isZeroLength() const {
            float sqlen = (m_data.x * m_data.x) + (m_data.y * m_data.y) + (m_data.z * m_data.z);
            return (sqlen < (Float_EPSILON * Float_EPSILON));
        }
        
        // 检查是否包含NaN值
        bool isNaN() const {
            return Math::isNan(m_data.x) || Math::isNan(m_data.y) || Math::isNan(m_data.z);
        }
        
        // 线性插值
        static Vector3 lerp(const Vector3& a, const Vector3& b, float t) {
            return Vector3(
                a.m_data.x + t * (b.m_data.x - a.m_data.x),
                a.m_data.y + t * (b.m_data.y - a.m_data.y),
                a.m_data.z + t * (b.m_data.z - a.m_data.z)
            );
        }
        
        // 静态常量
        static const Vector3 ZERO;
        static const Vector3 UNIT_X;
        static const Vector3 UNIT_Y;
        static const Vector3 UNIT_Z;
        static const Vector3 UNIT_SCALE;
        
        // 转换到GLM类型
        const glm::vec3& toGlm() const { return m_data; }
        glm::vec3& toGlm() { return m_data; }
    };
    
    // 静态常量定义
    inline const Vector3 Vector3::ZERO = Vector3(0.0f, 0.0f, 0.0f);
    inline const Vector3 Vector3::UNIT_X = Vector3(1.0f, 0.0f, 0.0f);
    inline const Vector3 Vector3::UNIT_Y = Vector3(0.0f, 1.0f, 0.0f);
    inline const Vector3 Vector3::UNIT_Z = Vector3(0.0f, 0.0f, 1.0f);
    inline const Vector3 Vector3::UNIT_SCALE = Vector3(1.0f, 1.0f, 1.0f);
}

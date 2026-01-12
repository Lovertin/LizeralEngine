#pragma once

#include "runtime/core/math/math.h"
#include "runtime/core/math/vector3.h"
#include "runtime/core/meta/reflection/reflection.h"
#include <glm/glm.hpp>
#include <cassert>

namespace Lizeral
{
    REFLECTION_TYPE(Vector4)
    CLASS(Vector4, Fields)
    {
        REFLECTION_BODY(Vector4);

    private:
        glm::vec4 m_data;

    public:
        Vector4() : m_data(0.0f, 0.0f, 0.0f, 0.0f) {}
        Vector4(float x_, float y_, float z_, float w_) : m_data(x_, y_, z_, w_) {}
        Vector4(const Vector3& v3, float w_) : m_data(v3.getX(), v3.getY(), v3.getZ(), w_) {}
        explicit Vector4(float coords[4]) : m_data(coords[0], coords[1], coords[2], coords[3]) {}
        explicit Vector4(const glm::vec4& v) : m_data(v) {}
        explicit Vector4(float scalar) : m_data(scalar, scalar, scalar, scalar) {}

        // 字段访问（用于反射）
        float getX() const { return m_data.x; }
        void setX(float value) { m_data.x = value; }
        
        float getY() const { return m_data.y; }
        void setY(float value) { m_data.y = value; }
        
        float getZ() const { return m_data.z; }
        void setZ(float value) { m_data.z = value; }
        
        float getW() const { return m_data.w; }
        void setW(float value) { m_data.w = value; }

        float operator[](size_t i) const
        {
            assert(i < 4);
            return m_data[i];
        }

        float& operator[](size_t i)
        {
            assert(i < 4);
            return m_data[i];
        }

        /// Pointer accessor for direct copying
        float* ptr() { return &m_data.x; }
        /// Pointer accessor for direct copying
        const float* ptr() const { return &m_data.x; }

        Vector4& operator=(float scalar)
        {
            m_data = glm::vec4(scalar, scalar, scalar, scalar);
            return *this;
        }

        bool operator==(const Vector4& rhs) const { return m_data == rhs.m_data; }
        bool operator!=(const Vector4& rhs) const { return m_data != rhs.m_data; }

        Vector4 operator+(const Vector4& rhs) const { return Vector4(m_data + rhs.m_data); }
        Vector4 operator-(const Vector4& rhs) const { return Vector4(m_data - rhs.m_data); }
        Vector4 operator*(float scalar) const { return Vector4(m_data * scalar); }
        Vector4 operator*(const Vector4& rhs) const { return Vector4(m_data * rhs.m_data); }
        Vector4 operator/(float scalar) const
        {
            assert(scalar != 0.0);
            return Vector4(m_data / scalar);
        }
        Vector4 operator/(const Vector4& rhs) const
        {
            assert(rhs.m_data.x != 0 && rhs.m_data.y != 0 && rhs.m_data.z != 0 && rhs.m_data.w != 0);
            return Vector4(m_data / rhs.m_data);
        }

        const Vector4& operator+() const { return *this; }
        Vector4 operator-() const { return Vector4(-m_data); }

        friend Vector4 operator*(float scalar, const Vector4& rhs)
        {
            return rhs * scalar;
        }

        friend Vector4 operator/(float scalar, const Vector4& rhs)
        {
            assert(rhs.m_data.x != 0 && rhs.m_data.y != 0 && rhs.m_data.z != 0 && rhs.m_data.w != 0);
            return Vector4(scalar / rhs.m_data.x, scalar / rhs.m_data.y, scalar / rhs.m_data.z, scalar / rhs.m_data.w);
        }

        friend Vector4 operator+(const Vector4& lhs, float rhs)
        {
            return Vector4(lhs.m_data + rhs);
        }

        friend Vector4 operator+(float lhs, const Vector4& rhs)
        {
            return rhs + lhs;
        }

        friend Vector4 operator-(const Vector4& lhs, float rhs)
        {
            return Vector4(lhs.m_data - rhs);
        }

        friend Vector4 operator-(float lhs, const Vector4& rhs)
        {
            return Vector4(lhs - rhs.m_data.x, lhs - rhs.m_data.y, lhs - rhs.m_data.z, lhs - rhs.m_data.w);
        }

        // arithmetic updates
        Vector4& operator+=(const Vector4& rhs)
        {
            m_data += rhs.m_data;
            return *this;
        }

        Vector4& operator-=(const Vector4& rhs)
        {
            m_data -= rhs.m_data;
            return *this;
        }

        Vector4& operator*=(float scalar)
        {
            m_data *= scalar;
            return *this;
        }

        Vector4& operator+=(float scalar)
        {
            m_data += scalar;
            return *this;
        }

        Vector4& operator-=(float scalar)
        {
            m_data -= scalar;
            return *this;
        }

        Vector4& operator*=(const Vector4& rhs)
        {
            m_data *= rhs.m_data;
            return *this;
        }

        Vector4& operator/=(float scalar)
        {
            assert(scalar != 0.0);
            m_data /= scalar;
            return *this;
        }

        Vector4& operator/=(const Vector4& rhs)
        {
            assert(rhs.m_data.x != 0 && rhs.m_data.y != 0 && rhs.m_data.z != 0 && rhs.m_data.w != 0);
            m_data /= rhs.m_data;
            return *this;
        }

        /** Calculates the dot (scalar) product of this vector with another.
        @param
        vec Vector with which to calculate the dot product (together
        with this one).
        @returns
        A float representing the dot product value.
        */
        float dotProduct(const Vector4& vec) const { return glm::dot(m_data, vec.m_data); }

        // 其他向量运算
        float length() const { return glm::length(m_data); }
        float squaredLength() const { return glm::dot(m_data, m_data); }
        
        void normalise() { m_data = glm::normalize(m_data); }
        Vector4 normalisedCopy() const { return Vector4(glm::normalize(m_data)); }
        
        /// Check whether this vector contains valid values
        bool isNaN() const { 
            return Math::isNan(m_data.x) || Math::isNan(m_data.y) || 
                   Math::isNan(m_data.z) || Math::isNan(m_data.w); 
        }

        // 线性插值
        static Vector4 lerp(const Vector4& a, const Vector4& b, float t) {
            return Vector4(
                a.m_data.x + t * (b.m_data.x - a.m_data.x),
                a.m_data.y + t * (b.m_data.y - a.m_data.y),
                a.m_data.z + t * (b.m_data.z - a.m_data.z),
                a.m_data.w + t * (b.m_data.w - a.m_data.w)
            );
        }

        // 转换为GLM类型
        const glm::vec4& toGlm() const { return m_data; }
        glm::vec4& toGlm() { return m_data; }

        // 转换为Vector3（忽略w分量）
        Vector3 xyz() const { return Vector3(m_data.x, m_data.y, m_data.z); }

        // special
        static const Vector4 ZERO;
        static const Vector4 UNIT_SCALE;
    };

    // 静态常量定义
    inline const Vector4 Vector4::ZERO = Vector4(0.0f, 0.0f, 0.0f, 0.0f);
    inline const Vector4 Vector4::UNIT_SCALE = Vector4(1.0f, 1.0f, 1.0f, 1.0f);
} // namespace Lizeral

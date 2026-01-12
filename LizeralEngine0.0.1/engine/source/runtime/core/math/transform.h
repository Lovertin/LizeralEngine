#pragma once
#include "matrix4.h"
#include "quaternion.h"
#include "vector3.h"
#include "../meta/reflection/reflection.h"

namespace Lizeral
{
    REFLECTION_TYPE(Transform)
    CLASS(Transform, Fields)
    {
        REFLECTION_BODY(Transform);

    public:
        Vector3    m_position {Vector3::ZERO};
        Vector3    m_scale {Vector3::UNIT_SCALE};
        Quaternion m_rotation {Quaternion::IDENTITY};

        Transform() = default;
        Transform(const Vector3& position, const Quaternion& rotation, const Vector3& scale) :
            m_position {position}, m_scale {scale}, m_rotation {rotation}
        {}

        Matrix4x4 getMatrix() const
        {
            Matrix4x4 temp;
            temp.makeTransform(m_position, m_scale, m_rotation);
            return temp;
        }

        // 获取逆矩阵
        Matrix4x4 getInverseMatrix() const
        {
            Matrix4x4 temp;
            temp.makeInverseTransform(m_position, m_scale, m_rotation);
            return temp;
        }

        // 组合变换
        Transform combine(const Transform& other) const
        {
            Transform result;
            // 先缩放，再旋转，最后平移
            result.m_scale = Vector3(
                m_scale.getX() * other.m_scale.getX(),
                m_scale.getY() * other.m_scale.getY(),
                m_scale.getZ() * other.m_scale.getZ()
            );
            
            // 旋转组合：先应用当前旋转，再应用其他旋转
            result.m_rotation = m_rotation * other.m_rotation;
            
            // 平移组合：先缩放旋转当前平移，再加上其他平移
            Vector3 rotatedScaledPosition = m_rotation * (m_scale * other.m_position);
            result.m_position = m_position + rotatedScaledPosition;
            
            return result;
        }

        // 逆变换
        Transform inverse() const
        {
            Transform result;
            // 逆缩放
            result.m_scale = Vector3(
                1.0f / m_scale.getX(),
                1.0f / m_scale.getY(),
                1.0f / m_scale.getZ()
            );
            
            // 逆旋转
            result.m_rotation = m_rotation.inverse();
            
            // 逆平移：先逆旋转，再逆缩放，最后取反
            Vector3 invPosition = -m_position;
            result.m_position = result.m_rotation * (result.m_scale * invPosition);
            
            return result;
        }

        // 变换点
        Vector3 transformPoint(const Vector3& point) const
        {
            // 应用缩放、旋转、平移
            Vector3 scaled = Vector3(
                point.getX() * m_scale.getX(),
                point.getY() * m_scale.getY(),
                point.getZ() * m_scale.getZ()
            );
            Vector3 rotated = m_rotation * scaled;
            return rotated + m_position;
        }

        // 变换方向（忽略平移）
        Vector3 transformDirection(const Vector3& direction) const
        {
            // 只应用缩放和旋转
            Vector3 scaled = Vector3(
                direction.getX() * m_scale.getX(),
                direction.getY() * m_scale.getY(),
                direction.getZ() * m_scale.getZ()
            );
            return m_rotation * scaled;
        }

        // 从矩阵设置变换
        void setFromMatrix(const Matrix4x4& matrix)
        {
            matrix.decomposition(m_position, m_scale, m_rotation);
        }

        // 重置为单位变换
        void reset()
        {
            m_position = Vector3::ZERO;
            m_scale = Vector3::UNIT_SCALE;
            m_rotation = Quaternion::IDENTITY;
        }

        // 检查是否为单位变换
        bool isIdentity() const
        {
            return m_position == Vector3::ZERO && 
                   m_scale == Vector3::UNIT_SCALE && 
                   m_rotation == Quaternion::IDENTITY;
        }

        // 插值
        static Transform lerp(const Transform& a, const Transform& b, float t)
        {
            Transform result;
            result.m_position = Vector3::lerp(a.m_position, b.m_position, t);
            result.m_scale = Vector3::lerp(a.m_scale, b.m_scale, t);
            result.m_rotation = Quaternion::nLerp(t, a.m_rotation, b.m_rotation);
            return result;
        }

        // 球面插值
        static Transform slerp(const Transform& a, const Transform& b, float t)
        {
            Transform result;
            result.m_position = Vector3::lerp(a.m_position, b.m_position, t);
            result.m_scale = Vector3::lerp(a.m_scale, b.m_scale, t);
            result.m_rotation = Quaternion::sLerp(t, a.m_rotation, b.m_rotation);
            return result;
        }

        // 运算符重载
        bool operator==(const Transform& other) const
        {
            return m_position == other.m_position && 
                   m_scale == other.m_scale && 
                   m_rotation == other.m_rotation;
        }

        bool operator!=(const Transform& other) const
        {
            return !(*this == other);
        }

        // 组合运算符
        Transform operator*(const Transform& other) const
        {
            return combine(other);
        }

        // 变换点运算符
        Vector3 operator*(const Vector3& point) const
        {
            return transformPoint(point);
        }
    };
} // namespace Lizeral

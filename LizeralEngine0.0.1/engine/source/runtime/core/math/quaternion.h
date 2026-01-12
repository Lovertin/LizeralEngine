#pragma once

#include "runtime/core/math/math.h"
#include "runtime/core/meta/reflection/reflection.h"
#include "runtime/core/math/vector3.h"  // 包含Vector3定义
#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/quaternion.hpp>
#include <cassert>

namespace Lizeral {
    // 前向声明
    class Matrix3x3;
    class Matrix4x4;

    REFLECTION_TYPE(Quaternion)
    CLASS(Quaternion, Fields) {
        REFLECTION_BODY(Quaternion);

    private:
        glm::quat m_data;  // GLM四元数：w + xi + yj + zk

    public:
        // 构造函数
        Quaternion() : m_data(1.0f, 0.0f, 0.0f, 0.0f) {}  // 单位四元数
        Quaternion(float w, float x, float y, float z) : m_data(w, x, y, z) {}
        
        // 从GLM四元数构造
        explicit Quaternion(const glm::quat& q) : m_data(q) {}
        
        // 从旋转矩阵构造（需要Matrix3x3实现）
        explicit Quaternion(const Matrix3x3& rot);
        
        // 从角度/轴构造
        Quaternion(const Radian& angle, const Vector3& axis);
        
        // 从3个正交轴构造
        Quaternion(const Vector3& xaxis, const Vector3& yaxis, const Vector3& zaxis);
        
        // 字段访问（用于反射）
        float getW() const { return m_data.w; }
        void setW(float value) { m_data.w = value; }
        
        float getX() const { return m_data.x; }
        void setX(float value) { m_data.x = value; }
        
        float getY() const { return m_data.y; }
        void setY(float value) { m_data.y = value; }
        
        float getZ() const { return m_data.z; }
        void setZ(float value) { m_data.z = value; }
        
        // 指针访问（用于直接复制）
        float* ptr() { return &m_data.x; }
        const float* ptr() const { return &m_data.x; }
        
        // 从旋转矩阵构造
        void fromRotationMatrix(const Matrix3x3& rotation);
        
        // 转换为旋转矩阵
        void toRotationMatrix(Matrix3x3& rotation) const;
        void toRotationMatrix(Matrix4x4& rotation) const;
        
        // 从角度/轴构造
        void fromAngleAxis(const Radian& angle, const Vector3& axis);
        
        // 静态方法：从角度/轴获取四元数
        static Quaternion getQuaternionFromAngleAxis(const Radian& angle, const Vector3& axis);
        
        // 从方向构造
        void fromDirection(const Vector3& direction, const Vector3& up_direction);
        
        // 静态方法：从方向获取四元数
        static Quaternion getQuaternionFromDirection(const Vector3& direction, const Vector3& up_direction);
        
        // 转换为角度/轴
        void toAngleAxis(Radian& angle, Vector3& axis) const;
        
        // 从3个正交轴构造
        void fromAxes(const Vector3& x_axis, const Vector3& y_axis, const Vector3& z_axis);
        
        // 获取3个正交轴
        void toAxes(Vector3& x_axis, Vector3& y_axis, Vector3& z_axis) const;
        
        // 获取局部轴
        Vector3 xAxis() const;
        Vector3 yAxis() const;
        Vector3 zAxis() const;
        
        // 基本运算符
        Quaternion operator+(const Quaternion& rhs) const {
            return Quaternion(m_data.w + rhs.m_data.w, 
                             m_data.x + rhs.m_data.x, 
                             m_data.y + rhs.m_data.y, 
                             m_data.z + rhs.m_data.z);
        }
        
        Quaternion operator-(const Quaternion& rhs) const {
            return Quaternion(m_data.w - rhs.m_data.w, 
                             m_data.x - rhs.m_data.x, 
                             m_data.y - rhs.m_data.y, 
                             m_data.z - rhs.m_data.z);
        }
        
        // 四元数乘法（使用GLM）
        Quaternion operator*(const Quaternion& rhs) const {
            return Quaternion(m_data * rhs.m_data);
        }
        
        Quaternion mul(const Quaternion& rhs) const { return (*this) * rhs; }
        
        // 标量乘法
        Quaternion operator*(float scalar) const {
            return Quaternion(m_data.w * scalar, 
                             m_data.x * scalar, 
                             m_data.y * scalar, 
                             m_data.z * scalar);
        }
        
        // 向量旋转（使用GLM）
        Vector3 operator*(const Vector3& rhs) const {
            glm::vec3 result = m_data * rhs.toGlm();
            return Vector3(result);
        }
        
        // 标量除法
        Quaternion operator/(float scalar) const {
            assert(scalar != 0.0f);
            return Quaternion(m_data.w / scalar, 
                             m_data.x / scalar, 
                             m_data.y / scalar, 
                             m_data.z / scalar);
        }
        
        // 从角度/轴构造的实现
        void fromAngleAxis(const Radian& angle, const Vector3& axis) {
            float halfAngle = angle.valueRadians() * 0.5f;
            float s = std::sin(halfAngle);
            m_data.w = std::cos(halfAngle);
            m_data.x = axis.getX() * s;
            m_data.y = axis.getY() * s;
            m_data.z = axis.getZ() * s;
        }
        
        // 转换为角度/轴
        void toAngleAxis(Radian& angle, Vector3& axis) const {
            // 检查四元数是否为单位四元数
            if (m_data.w > 1.0f) {
                Quaternion q = *this;
                q.normalise();
                q.toAngleAxis(angle, axis);
                return;
            }
            
            float squaredSin = 1.0f - m_data.w * m_data.w;
            if (squaredSin > 0.0f) {
                float s = std::sqrt(squaredSin);
                angle = Radian(2.0f * std::acos(m_data.w));
                axis.setX(m_data.x / s);
                axis.setY(m_data.y / s);
                axis.setZ(m_data.z / s);
            } else {
                // 对于单位四元数，角度为0，轴任意
                angle = Radian(0.0f);
                axis.setX(1.0f);
                axis.setY(0.0f);
                axis.setZ(0.0f);
            }
        }
        
        // 获取局部轴
        Vector3 xAxis() const {
            float fTy  = 2.0f * m_data.y;
            float fTz  = 2.0f * m_data.z;
            float fTwy = fTy * m_data.w;
            float fTwz = fTz * m_data.w;
            float fTxy = fTy * m_data.x;
            float fTxz = fTz * m_data.x;
            float fTyy = fTy * m_data.y;
            float fTzz = fTz * m_data.z;
            
            return Vector3(1.0f - (fTyy + fTzz), fTxy + fTwz, fTxz - fTwy);
        }
        
        Vector3 yAxis() const {
            float fTx  = 2.0f * m_data.x;
            float fTy  = 2.0f * m_data.y;
            float fTz  = 2.0f * m_data.z;
            float fTwx = fTx * m_data.w;
            float fTwz = fTz * m_data.w;
            float fTxx = fTx * m_data.x;
            float fTxy = fTy * m_data.x;
            float fTyz = fTz * m_data.y;
            float fTzz = fTz * m_data.z;
            
            return Vector3(fTxy - fTwz, 1.0f - (fTxx + fTzz), fTyz + fTwx);
        }
        
        Vector3 zAxis() const {
            float fTx  = 2.0f * m_data.x;
            float fTy  = 2.0f * m_data.y;
            float fTz  = 2.0f * m_data.z;
            float fTwx = fTx * m_data.w;
            float fTwy = fTy * m_data.w;
            float fTxx = fTx * m_data.x;
            float fTxz = fTz * m_data.x;
            float fTyy = fTy * m_data.y;
            float fTyz = fTz * m_data.y;
            
            return Vector3(fTxz + fTwy, fTyz - fTwx, 1.0f - (fTxx + fTyy));
        }
        
        // 从方向构造
        void fromDirection(const Vector3& direction, const Vector3& up_direction) {
            Vector3 forward = direction.normalisedCopy();
            Vector3 up = up_direction.normalisedCopy();
            Vector3 right = forward.cross(up).normalisedCopy();
            up = right.cross(forward);
            
            fromAxes(right, up, -forward);
        }
        
        // 从3个正交轴构造
        void fromAxes(const Vector3& x_axis, const Vector3& y_axis, const Vector3& z_axis) {
            Matrix3x3 rot;
            // 这里需要Matrix3x3的实现
            // 暂时留空，等Matrix3x3实现后再完成
        }
        
        // 获取3个正交轴
        void toAxes(Vector3& x_axis, Vector3& y_axis, Vector3& z_axis) const {
            x_axis = xAxis();
            y_axis = yAxis();
            z_axis = zAxis();
        }
        
        // 静态方法：从角度/轴获取四元数
        static Quaternion getQuaternionFromAngleAxis(const Radian& angle, const Vector3& axis) {
            Quaternion q;
            q.fromAngleAxis(angle, axis);
            return q;
        }
        
        // 静态方法：从方向获取四元数
        static Quaternion getQuaternionFromDirection(const Vector3& direction, const Vector3& up_direction) {
            Quaternion q;
            q.fromDirection(direction, up_direction);
            return q;
        }
        
        // 友元函数：标量乘法
        friend Quaternion operator*(float scalar, const Quaternion& rhs) {
            return rhs * scalar;
        }
        
        // 取反
        Quaternion operator-() const {
            return Quaternion(-m_data.w, -m_data.x, -m_data.y, -m_data.z);
        }
        
        // 相等比较
        bool operator==(const Quaternion& rhs) const {
            return (m_data.x == rhs.m_data.x) && 
                   (m_data.y == rhs.m_data.y) && 
                   (m_data.z == rhs.m_data.z) && 
                   (m_data.w == rhs.m_data.w);
        }
        
        bool operator!=(const Quaternion& rhs) const {
            return !(*this == rhs);
        }
        
        // 检查是否包含NaN值
        bool isNaN() const {
            return Math::isNan(m_data.x) || Math::isNan(m_data.y) || 
                   Math::isNan(m_data.z) || Math::isNan(m_data.w);
        }
        
        // 点积
        float dot(const Quaternion& rkQ) const {
            return glm::dot(m_data, rkQ.m_data);
        }
        
        // 长度
        float length() const {
            return glm::length(m_data);
        }
        
        // 归一化
        void normalise() {
            m_data = glm::normalize(m_data);
        }
        
        // 逆四元数
        Quaternion inverse() const {
            return Quaternion(glm::inverse(m_data));
        }
        
        // 共轭四元数
        Quaternion conjugate() const {
            return Quaternion(glm::conjugate(m_data));
        }
        
        // 获取欧拉角（需要实现）
        Radian getRoll(bool reproject_axis = true) const;
        Radian getPitch(bool reproject_axis = true) const;
        Radian getYaw(bool reproject_axis = true) const;
        
        // 球面线性插值（使用GLM）
        static Quaternion sLerp(float t, const Quaternion& kp, const Quaternion& kq, bool shortest_path = false) {
            if (shortest_path) {
                // 使用GLM的slerp，确保最短路径
                return Quaternion(glm::slerp(kp.m_data, kq.m_data, t));
            } else {
                // 标准slerp
                return Quaternion(glm::mix(kp.m_data, kq.m_data, t));
            }
        }
        
        // 归一化线性插值（使用GLM）
        static Quaternion nLerp(float t, const Quaternion& kp, const Quaternion& kq, bool shortest_path = false) {
            // GLM的mix函数已经做了归一化线性插值
            glm::quat result = glm::mix(kp.m_data, kq.m_data, t);
            result = glm::normalize(result);
            return Quaternion(result);
        }
        
        // 转换到GLM类型
        const glm::quat& toGlm() const { return m_data; }
        glm::quat& toGlm() { return m_data; }
        
        // 静态常量
        static const Quaternion ZERO;
        static const Quaternion IDENTITY;
        
        static const float k_epsilon;
    };
    
    // 静态常量定义
    inline const Quaternion Quaternion::ZERO = Quaternion(0.0f, 0.0f, 0.0f, 0.0f);
    inline const Quaternion Quaternion::IDENTITY = Quaternion(1.0f, 0.0f, 0.0f, 0.0f);
    inline const float Quaternion::k_epsilon = 1e-6f;
}

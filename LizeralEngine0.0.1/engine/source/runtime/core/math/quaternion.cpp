#include "runtime/core/math/quaternion.h"
#include "runtime/core/math/matrix3.h"
#include "runtime/core/math/matrix4.h"
#include "runtime/core/math/vector3.h"

namespace Lizeral
{
    // 构造函数实现
    Quaternion::Quaternion(const Matrix3x3& rot)
    {
        fromRotationMatrix(rot);
    }
    
    Quaternion::Quaternion(const Radian& angle, const Vector3& axis)
    {
        fromAngleAxis(angle, axis);
    }
    
    Quaternion::Quaternion(const Vector3& xaxis, const Vector3& yaxis, const Vector3& zaxis)
    {
        fromAxes(xaxis, yaxis, zaxis);
    }
    
    //-----------------------------------------------------------------------
    void Quaternion::fromRotationMatrix(const Matrix3x3& rotation)
    {
        // Algorithm in Ken Shoemake's article in 1987 SIGGRAPH course notes
        // article "Quaternion Calculus and Fast Animation".

        float trace = rotation[0][0] + rotation[1][1] + rotation[2][2];
        float root;

        if (trace > 0.0f)
        {
            // |w| > 1/2, may as well choose w > 1/2
            root = std::sqrt(trace + 1.0f); // 2w
            m_data.w    = 0.5f * root;
            root = 0.5f / root; // 1/(4w)
            m_data.x    = (rotation[2][1] - rotation[1][2]) * root;
            m_data.y    = (rotation[0][2] - rotation[2][0]) * root;
            m_data.z    = (rotation[1][0] - rotation[0][1]) * root;
        }
        else
        {
            // |w| <= 1/2
            size_t s_iNext[3] = {1, 2, 0};
            size_t i          = 0;
            if (rotation[1][1] > rotation[0][0])
                i = 1;
            if (rotation[2][2] > rotation[i][i])
                i = 2;
            size_t j = s_iNext[i];
            size_t k = s_iNext[j];

            root              = std::sqrt(rotation[i][i] - rotation[j][j] - rotation[k][k] + 1.0f);
            float* apkQuat[3] = {&m_data.x, &m_data.y, &m_data.z};
            *apkQuat[i]       = 0.5f * root;
            root              = 0.5f / root;
            m_data.w                 = (rotation[k][j] - rotation[j][k]) * root;
            *apkQuat[j]       = (rotation[j][i] + rotation[i][j]) * root;
            *apkQuat[k]       = (rotation[k][i] + rotation[i][k]) * root;
        }
    }
    
    //-----------------------------------------------------------------------
    void Quaternion::toRotationMatrix(Matrix3x3& kRot) const
    {
        float fTx  = m_data.x + m_data.x;   // 2x
        float fTy  = m_data.y + m_data.y;   // 2y
        float fTz  = m_data.z + m_data.z;   // 2z
        float fTwx = fTx * m_data.w; // 2xw
        float fTwy = fTy * m_data.w; // 2yw
        float fTwz = fTz * m_data.w; // 2z2
        float fTxx = fTx * m_data.x; // 2x^2
        float fTxy = fTy * m_data.x; // 2xy
        float fTxz = fTz * m_data.x; // 2xz
        float fTyy = fTy * m_data.y; // 2y^2
        float fTyz = fTz * m_data.y; // 2yz
        float fTzz = fTz * m_data.z; // 2z^2

        kRot[0][0] = 1.0f - (fTyy + fTzz); // 1 - 2y^2 - 2z^2
        kRot[0][1] = fTxy - fTwz;          // 2xy - 2wz
        kRot[0][2] = fTxz + fTwy;          // 2xz + 2wy
        kRot[1][0] = fTxy + fTwz;          // 2xy + 2wz
        kRot[1][1] = 1.0f - (fTxx + fTzz); // 1 - 2x^2 - 2z^2
        kRot[1][2] = fTyz - fTwx;          // 2yz - 2wx
        kRot[2][0] = fTxz - fTwy;          // 2xz - 2wy
        kRot[2][1] = fTyz + fTwx;          // 2yz + 2wx
        kRot[2][2] = 1.0f - (fTxx + fTyy); // 1 - 2x^2 - 2y^2
    }

    void Quaternion::toRotationMatrix(Matrix4x4& kRot) const
    {
        float fTx  = m_data.x + m_data.x;   // 2x
        float fTy  = m_data.y + m_data.y;   // 2y
        float fTz  = m_data.z + m_data.z;   // 2z
        float fTwx = fTx * m_data.w; // 2xw
        float fTwy = fTy * m_data.w; // 2yw
        float fTwz = fTz * m_data.w; // 2z2
        float fTxx = fTx * m_data.x; // 2x^2
        float fTxy = fTy * m_data.x; // 2xy
        float fTxz = fTz * m_data.x; // 2xz
        float fTyy = fTy * m_data.y; // 2y^2
        float fTyz = fTz * m_data.y; // 2yz
        float fTzz = fTz * m_data.z; // 2z^2

        kRot[0][0] = 1.0f - (fTyy + fTzz); // 1 - 2y^2 - 2z^2
        kRot[0][1] = fTxy - fTwz;          // 2xy - 2wz
        kRot[0][2] = fTxz + fTwy;          // 2xz + 2wy
        kRot[0][3] = 0;
        kRot[1][0] = fTxy + fTwz;          // 2xy + 2wz
        kRot[1][1] = 1.0f - (fTxx + fTzz); // 1 - 2x^2 - 2z^2
        kRot[1][2] = fTyz - fTwx;          // 2yz - 2wx
        kRot[1][3] = 0;
        kRot[2][0] = fTxz - fTwy;          // 2xz - 2wy
        kRot[2][1] = fTyz + fTwx;          // 2yz + 2wx
        kRot[2][2] = 1.0f - (fTxx + fTyy); // 1 - 2x^2 - 2y^2
        kRot[2][3] = 0;
        kRot[3][0] = 0;
        kRot[3][1] = 0;
        kRot[3][2] = 0;
        kRot[3][3] = 1;
    }

    //-----------------------------------------------------------------------
    void Quaternion::fromAngleAxis(const Radian& angle, const Vector3& axis)
    {
        // ASSERT:  axis[] is unit length
        //
        // The quaternion representing the rotation is
        //   q = cos(A/2)+sin(A/2)*(x*i+y*j+z*k)
        Radian half_angle(0.5f * angle.valueRadians());
        float  sin_v = Math::sin(half_angle);
        m_data.w            = Math::cos(half_angle);
        m_data.x            = sin_v * axis.getX();
        m_data.y            = sin_v * axis.getY();
        m_data.z            = sin_v * axis.getZ();
    }

    Quaternion Quaternion::getQuaternionFromAngleAxis(const Radian& angle, const Vector3& axis)
    {
        Quaternion q;
        q.fromAngleAxis(angle, axis);
        return q;
    }

    void Quaternion::fromDirection(const Vector3& direction, const Vector3& up_direction)
    {
        Vector3 forward_direction = direction;
        forward_direction.setZ(0.0f);
        forward_direction.normalise();

        Vector3 left_direction = up_direction.cross(forward_direction);

        fromAxes(left_direction, -forward_direction, up_direction);
        normalise();
    }

    Quaternion Quaternion::getQuaternionFromDirection(const Vector3& direction, const Vector3& up_direction)
    {
        Quaternion object_orientation;
        object_orientation.fromDirection(direction, up_direction);
        return object_orientation;
    }

    void Quaternion::toAngleAxis(Radian& angle, Vector3& axis) const
    {
        // The quaternion representing the rotation is
        //   q = cos(A/2)+sin(A/2)*(x*i+y*j+z*k)

        float sqr_len = m_data.x * m_data.x + m_data.y * m_data.y + m_data.z * m_data.z;
        if (sqr_len > 0.0f)
        {
            angle         = Radian(2.0f * Math::acos(m_data.w));
            float inv_len = Math::invSqrt(sqr_len);
            axis.setX(m_data.x * inv_len);
            axis.setY(m_data.y * inv_len);
            axis.setZ(m_data.z * inv_len);
        }
        else
        {
            // angle is 0 (mod 2*pi), so any axis will do
            angle  = Radian(0.0f);
            axis.setX(1.0f);
            axis.setY(0.0f);
            axis.setZ(0.0f);
        }
    }

    //-----------------------------------------------------------------------
    void Quaternion::fromAxes(const Vector3& xaxis, const Vector3& yaxis, const Vector3& zaxis)
    {
        Matrix3x3 rot;

        rot[0][0] = xaxis.getX();
        rot[1][0] = xaxis.getY();
        rot[2][0] = xaxis.getZ();

        rot[0][1] = yaxis.getX();
        rot[1][1] = yaxis.getY();
        rot[2][1] = yaxis.getZ();

        rot[0][2] = zaxis.getX();
        rot[1][2] = zaxis.getY();
        rot[2][2] = zaxis.getZ();

        fromRotationMatrix(rot);
    }
    
    //-----------------------------------------------------------------------
    Vector3 Quaternion::xAxis() const
    {
        // float tx  = 2.0*x;
        float ty  = 2.0f * m_data.y;
        float tz  = 2.0f * m_data.z;
        float twy = ty * m_data.w;
        float twz = tz * m_data.w;
        float txy = ty * m_data.x;
        float txz = tz * m_data.x;
        float tyy = ty * m_data.y;
        float tzz = tz * m_data.z;

        return Vector3(1.0f - (tyy + tzz), txy + twz, txz - twy);
    }
    
    //-----------------------------------------------------------------------
    Vector3 Quaternion::yAxis() const
    {
        float tx  = 2.0f * m_data.x;
        float ty  = 2.0f * m_data.y;
        float tz  = 2.0f * m_data.z;
        float twx = tx * m_data.w;
        float twz = tz * m_data.w;
        float txx = tx * m_data.x;
        float txy = ty * m_data.x;
        float tyz = tz * m_data.y;
        float tzz = tz * m_data.z;

        return Vector3(txy - twz, 1.0f - (txx + tzz), tyz + twx);
    }
    
    //-----------------------------------------------------------------------
    Vector3 Quaternion::zAxis() const
    {
        float tx  = 2.0f * m_data.x;
        float ty  = 2.0f * m_data.y;
        float tz  = 2.0f * m_data.z;
        float twx = tx * m_data.w;
        float twy = ty * m_data.w;
        float txx = tx * m_data.x;
        float txz = tz * m_data.x;
        float tyy = ty * m_data.y;
        float tyz = tz * m_data.y;

        return Vector3(txz + twy, tyz - twx, 1.0f - (txx + tyy));
    }
    
    //-----------------------------------------------------------------------
    void Quaternion::toAxes(Vector3& xaxis, Vector3& yaxis, Vector3& zaxis) const
    {
        Matrix3x3 rot;

        toRotationMatrix(rot);

        xaxis.setX(rot[0][0]);
        xaxis.setY(rot[1][0]);
        xaxis.setZ(rot[2][0]);

        yaxis.setX(rot[0][1]);
        yaxis.setY(rot[1][1]);
        yaxis.setZ(rot[2][1]);

        zaxis.setX(rot[0][2]);
        zaxis.setY(rot[1][2]);
        zaxis.setZ(rot[2][2]);
    }

    Vector3 Quaternion::operator*(const Vector3& v) const
    {
        // nVidia SDK implementation
        Vector3 uv, uuv;
        Vector3 qvec(m_data.x, m_data.y, m_data.z);
        uv  = qvec.cross(v);
        uuv = qvec.cross(uv);
        uv *= (2.0f * m_data.w);
        uuv *= 2.0f;

        return v + uv + uuv;
    }

    Radian Quaternion::getYaw(bool reproject_axis) const
    {
        if (reproject_axis)
        {
            // roll = atan2(localx.y, localx.x)
            // pick parts of xAxis() implementation that we need
            //  float tx  = 2.0*x;
            float ty  = 2.0f * m_data.y;
            float tz  = 2.0f * m_data.z;
            float twz = tz * m_data.w;
            float txy = ty * m_data.x;
            float tyy = ty * m_data.y;
            float tzz = tz * m_data.z;

            return Radian(Math::atan2(txy + twz, 1.0f - (tyy + tzz)));
        }
        else
        {
            return Radian(Math::atan2(2 * (m_data.x * m_data.y + m_data.w * m_data.z), 
                                      m_data.w * m_data.w + m_data.x * m_data.x - m_data.y * m_data.y - m_data.z * m_data.z));
        }
    }
    
    //-----------------------------------------------------------------------
    Radian Quaternion::getPitch(bool reproject_axis) const
    {
        if (reproject_axis)
        {
            // pitch = atan2(localy.z, localy.y)
            // pick parts of yAxis() implementation that we need
            float tx = 2.0f * m_data.x;
            //  float ty  = 2.0f*y;
            float tz  = 2.0f * m_data.z;
            float twx = tx * m_data.w;
            float txx = tx * m_data.x;
            float tyz = tz * m_data.y;
            float tzz = tz * m_data.z;

            return Radian(Math::atan2(tyz + twx, 1.0f - (txx + tzz)));
        }
        else
        {
            // internal version
            return Radian(Math::atan2(2 * (m_data.y * m_data.z + m_data.w * m_data.x), 
                                      m_data.w * m_data.w - m_data.x * m_data.x - m_data.y * m_data.y + m_data.z * m_data.z));
        }
    }
    
    //-----------------------------------------------------------------------
    Radian Quaternion::getRoll(bool reproject_axis) const
    {
        if (reproject_axis)
        {
            // yaw = atan2(localz.x, localz.z)
            // pick parts of zAxis() implementation that we need
            float tx  = 2.0f * m_data.x;
            float ty  = 2.0f * m_data.y;
            float tz  = 2.0f * m_data.z;
            float twy = ty * m_data.w;
            float txx = tx * m_data.x;
            float txz = tz * m_data.x;
            float tyy = ty * m_data.y;

            return Radian(Math::atan2(txz + twy, 1.0f - (txx + tyy)));
        }
        else
        {
            // internal version
            return Radian(Math::asin(-2 * (m_data.x * m_data.z - m_data.w * m_data.y)));
        }
    }

    //-----------------------------------------------------------------------
    Quaternion Quaternion::sLerp(float t, const Quaternion& kp, const Quaternion& kq, bool shortest_path)
    {
        float      cos_v = kp.dot(kq);
        Quaternion kt;

        // Do we need to invert rotation?
        if (cos_v < 0.0f && shortest_path)
        {
            cos_v = -cos_v;
            kt    = -kq;
        }
        else
        {
            kt = kq;
        }

        if (Math::abs(cos_v) < 1 - k_epsilon)
        {
            // Standard case (slerp)
            float  sin_v   = Math::sqrt(1 - Math::sqr(cos_v));
            Radian angle   = Math::atan2(sin_v, cos_v);
            float  inv_sin = 1.0f / sin_v;
            float  coeff0  = Math::sin((1.0f - t) * angle.valueRadians()) * inv_sin;
            float  coeff1  = Math::sin(t * angle.valueRadians()) * inv_sin;
            return coeff0 * kp + coeff1 * kt;
        }
        else
        {
            // There are two situations:
            // 1. "rkP" and "rkQ" are very close (fCos ~= +1), so we can do a linear
            //    interpolation safely.
            // 2. "rkP" and "rkQ" are almost inverse of each other (fCos ~= -1), there
            //    are an infinite number of possibilities interpolation. but we haven't
            //    have method to fix this case, so just use linear interpolation here.
            Quaternion r = (1.0f - t) * kp + t * kt;
            // taking the complement requires renormalization
            r.normalise();
            return r;
        }
    }

    //-----------------------------------------------------------------------
    Quaternion Quaternion::nLerp(float t, const Quaternion& kp, const Quaternion& kq, bool shortest_path)
    {
        Quaternion result;
        float      cos_value = kp.dot(kq);
        if (cos_value < 0.0f && shortest_path)
        {
            result = kp + t * ((-kq) - kp);
        }
        else
        {
            result = kp + t * (kq - kp);
        }
        result.normalise();
        return result;
    }
} // namespace Lizeral

#include "runtime/core/math/matrix3.h"
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/matrix_decompose.hpp>

namespace Lizeral
{
    const Matrix3x3 Matrix3x3::ZERO(0, 0, 0, 0, 0, 0, 0, 0, 0);
    const Matrix3x3 Matrix3x3::IDENTITY(1, 0, 0, 0, 1, 0, 0, 0, 1);

    //-----------------------------------------------------------------------
    void Matrix3x3::setColumn(size_t col_index, const Vector3& vec)
    {
        m_mat[0][col_index] = vec.getX();
        m_mat[1][col_index] = vec.getY();
        m_mat[2][col_index] = vec.getZ();
    }
    
    //-----------------------------------------------------------------------
    void Matrix3x3::fromAxes(const Vector3& x_axis, const Vector3& y_axis, const Vector3& z_axis)
    {
        setColumn(0, x_axis);
        setColumn(1, y_axis);
        setColumn(2, z_axis);
    }

    void Matrix3x3::calculateQDUDecomposition(Matrix3x3& out_Q, Vector3& out_D, Vector3& out_U) const
    {
        // 使用GLM的QR分解
        glm::mat3 glm_mat(
            m_mat[0][0], m_mat[1][0], m_mat[2][0],
            m_mat[0][1], m_mat[1][1], m_mat[2][1],
            m_mat[0][2], m_mat[1][2], m_mat[2][2]
        );
        
        // GLM没有直接的QR分解，我们使用原始算法
        // 但为了简化，我们可以使用特征值分解或SVD
        // 这里暂时使用原始算法
        
        // Factor M = QR = QDU where Q is orthogonal, D is diagonal,
        // and U is upper triangular with ones on its diagonal.
        
        // build orthogonal matrix Q
        float inv_length = m_mat[0][0] * m_mat[0][0] + m_mat[1][0] * m_mat[1][0] + m_mat[2][0] * m_mat[2][0];
        if (!Math::realEqual(inv_length, 0))
            inv_length = Math::invSqrt(inv_length);

        out_Q[0][0] = m_mat[0][0] * inv_length;
        out_Q[1][0] = m_mat[1][0] * inv_length;
        out_Q[2][0] = m_mat[2][0] * inv_length;

        float dot   = out_Q[0][0] * m_mat[0][1] + out_Q[1][0] * m_mat[1][1] + out_Q[2][0] * m_mat[2][1];
        out_Q[0][1] = m_mat[0][1] - dot * out_Q[0][0];
        out_Q[1][1] = m_mat[1][1] - dot * out_Q[1][0];
        out_Q[2][1] = m_mat[2][1] - dot * out_Q[2][0];
        inv_length  = out_Q[0][1] * out_Q[0][1] + out_Q[1][1] * out_Q[1][1] + out_Q[2][1] * out_Q[2][1];
        if (!Math::realEqual(inv_length, 0))
            inv_length = Math::invSqrt(inv_length);

        out_Q[0][1] *= inv_length;
        out_Q[1][1] *= inv_length;
        out_Q[2][1] *= inv_length;

        dot         = out_Q[0][0] * m_mat[0][2] + out_Q[1][0] * m_mat[1][2] + out_Q[2][0] * m_mat[2][2];
        out_Q[0][2] = m_mat[0][2] - dot * out_Q[0][0];
        out_Q[1][2] = m_mat[1][2] - dot * out_Q[1][0];
        out_Q[2][2] = m_mat[2][2] - dot * out_Q[2][0];
        dot         = out_Q[0][1] * m_mat[0][2] + out_Q[1][1] * m_mat[1][2] + out_Q[2][1] * m_mat[2][2];
        out_Q[0][2] -= dot * out_Q[0][1];
        out_Q[1][2] -= dot * out_Q[1][1];
        out_Q[2][2] -= dot * out_Q[2][1];
        inv_length = out_Q[0][2] * out_Q[0][2] + out_Q[1][2] * out_Q[1][2] + out_Q[2][2] * out_Q[2][2];
        if (!Math::realEqual(inv_length, 0))
            inv_length = Math::invSqrt(inv_length);

        out_Q[0][2] *= inv_length;
        out_Q[1][2] *= inv_length;
        out_Q[2][2] *= inv_length;

        // guarantee that orthogonal matrix has determinant 1 (no reflections)
        float det = out_Q[0][0] * out_Q[1][1] * out_Q[2][2] + out_Q[0][1] * out_Q[1][2] * out_Q[2][0] +
                    out_Q[0][2] * out_Q[1][0] * out_Q[2][1] - out_Q[0][2] * out_Q[1][1] * out_Q[2][0] -
                    out_Q[0][1] * out_Q[1][0] * out_Q[2][2] - out_Q[0][0] * out_Q[1][2] * out_Q[2][1];

        if (det < 0.0)
        {
            for (size_t row_index = 0; row_index < 3; row_index++)
                for (size_t rol_index = 0; rol_index < 3; rol_index++)
                    out_Q[row_index][rol_index] = -out_Q[row_index][rol_index];
        }

        // build "right" matrix R
        Matrix3x3 R;
        R[0][0] = out_Q[0][0] * m_mat[0][0] + out_Q[1][0] * m_mat[1][0] + out_Q[2][0] * m_mat[2][0];
        R[0][1] = out_Q[0][0] * m_mat[0][1] + out_Q[1][0] * m_mat[1][1] + out_Q[2][0] * m_mat[2][1];
        R[1][1] = out_Q[0][1] * m_mat[0][1] + out_Q[1][1] * m_mat[1][1] + out_Q[2][1] * m_mat[2][1];
        R[0][2] = out_Q[0][0] * m_mat[0][2] + out_Q[1][0] * m_mat[1][2] + out_Q[2][0] * m_mat[2][2];
        R[1][2] = out_Q[0][1] * m_mat[0][2] + out_Q[1][1] * m_mat[1][2] + out_Q[2][1] * m_mat[2][2];
        R[2][2] = out_Q[0][2] * m_mat[0][2] + out_Q[1][2] * m_mat[1][2] + out_Q[2][2] * m_mat[2][2];

        // the scaling component
        out_D[0] = R[0][0];
        out_D[1] = R[1][1];
        out_D[2] = R[2][2];

        // the shear component
        float inv_d0 = 1.0f / out_D[0];
        out_U[0]     = R[0][1] * inv_d0;
        out_U[1]     = R[0][2] * inv_d0;
        out_U[2]     = R[1][2] / out_D[1];
    }

    void Matrix3x3::toAngleAxis(Vector3& axis, Radian& radian) const
    {
        // 使用GLM将旋转矩阵转换为角度/轴
        glm::mat3 glm_mat(
            m_mat[0][0], m_mat[1][0], m_mat[2][0],
            m_mat[0][1], m_mat[1][1], m_mat[2][1],
            m_mat[0][2], m_mat[1][2], m_mat[2][2]
        );
        
        // GLM没有直接的矩阵到角度/轴的转换，我们使用原始算法
        float trace = m_mat[0][0] + m_mat[1][1] + m_mat[2][2];
        float cos_v = 0.5f * (trace - 1.0f);
        radian      = Math::acos(cos_v); // in [0,PI]

        if (radian > Radian(0.0))
        {
            if (radian < Radian(Math_PI))
            {
                axis.setX(m_mat[2][1] - m_mat[1][2]);
                axis.setY(m_mat[0][2] - m_mat[2][0]);
                axis.setZ(m_mat[1][0] - m_mat[0][1]);
                axis.normalise();
            }
            else
            {
                // angle is PI
                float half_inv;
                if (m_mat[0][0] >= m_mat[1][1])
                {
                    // r00 >= r11
                    if (m_mat[0][0] >= m_mat[2][2])
                    {
                        // r00 is maximum diagonal term
                        axis.setX(0.5f * Math::sqrt(m_mat[0][0] - m_mat[1][1] - m_mat[2][2] + 1.0f));
                        half_inv = 0.5f / axis.getX();
                        axis.setY(half_inv * m_mat[0][1]);
                        axis.setZ(half_inv * m_mat[0][2]);
                    }
                    else
                    {
                        // r22 is maximum diagonal term
                        axis.setZ(0.5f * Math::sqrt(m_mat[2][2] - m_mat[0][0] - m_mat[1][1] + 1.0f));
                        half_inv = 0.5f / axis.getZ();
                        axis.setX(half_inv * m_mat[0][2]);
                        axis.setY(half_inv * m_mat[1][2]);
                    }
                }
                else
                {
                    // r11 > r00
                    if (m_mat[1][1] >= m_mat[2][2])
                    {
                        // r11 is maximum diagonal term
                        axis.setY(0.5f * Math::sqrt(m_mat[1][1] - m_mat[0][0] - m_mat[2][2] + 1.0f));
                        half_inv = 0.5f / axis.getY();
                        axis.setX(half_inv * m_mat[0][1]);
                        axis.setZ(half_inv * m_mat[1][2]);
                    }
                    else
                    {
                        // r22 is maximum diagonal term
                        axis.setZ(0.5f * Math::sqrt(m_mat[2][2] - m_mat[0][0] - m_mat[1][1] + 1.0f));
                        half_inv = 0.5f / axis.getZ();
                        axis.setX(half_inv * m_mat[0][2]);
                        axis.setY(half_inv * m_mat[1][2]);
                    }
                }
            }
        }
        else
        {
            // The angle is 0 and the matrix is the identity.  Any axis will
            // work, so just use the x-axis.
            axis.setX(1.0);
            axis.setY(0.0);
            axis.setZ(0.0);
        }
    }
    
    //-----------------------------------------------------------------------
    void Matrix3x3::fromAngleAxis(const Vector3& axis, const Radian& radian)
    {
        // 使用GLM从角度/轴创建旋转矩阵
        // 对于3x3矩阵，我们可以使用四元数转换
        glm::vec3 glm_axis(axis.getX(), axis.getY(), axis.getZ());
        glm::quat glm_quat = glm::angleAxis(radian.valueRadians(), glm_axis);
        glm::mat3 glm_mat = glm::mat3_cast(glm_quat);
        
        // 将GLM矩阵复制到我们的矩阵
        m_mat[0][0] = glm_mat[0][0]; m_mat[0][1] = glm_mat[1][0]; m_mat[0][2] = glm_mat[2][0];
        m_mat[1][0] = glm_mat[0][1]; m_mat[1][1] = glm_mat[1][1]; m_mat[1][2] = glm_mat[2][1];
        m_mat[2][0] = glm_mat[0][2]; m_mat[2][1] = glm_mat[1][2]; m_mat[2][2] = glm_mat[2][2];
    }
} // namespace Lizeral

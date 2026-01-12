#pragma once

#include "runtime/core/math/math.h"
#include "runtime/core/math/matrix3.h"
#include "runtime/core/math/quaternion.h"
#include "runtime/core/math/vector3.h"
#include "runtime/core/math/vector4.h"
#include "runtime/core/meta/reflection/reflection.h"
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/matrix_decompose.hpp>
#include <cassert>

namespace Lizeral
{
    /** Class encapsulating a standard 4x4 homogeneous matrix.
    @remarks
    Lizeral uses column vectors when applying matrix multiplications,
    This means a vector is represented as a single column, 4-row
    matrix. This has the effect that the transformations implemented
    by the matrices happens right-to-left e.g. if vector V is to be
    transformed by M1 then M2 then M3, the calculation would be
    M3 * M2 * M1 * V. The order that matrices are concatenated is
    vital since matrix multiplication is not commutative, i.e. you
    can get a different result if you concatenate in the wrong order.
    @par
    The use of column vectors and right-to-left ordering is the
    standard in most mathematical texts, and is the same as used in
    OpenGL.
    */
    
    // 反射类型，用于序列化
    REFLECTION_TYPE(Matrix4x4_)
    CLASS(Matrix4x4_, Fields)
    {
        REFLECTION_BODY(Matrix4x4_);

    public:
        Matrix4x4_() {}
        float v0 {1.f};
        float v1 {0};
        float v2 {0};
        float v3 {0};
        float v4 {0};
        float v5 {1.f};
        float v6 {0};
        float v7 {0};
        float v8 {0};
        float v9 {0};
        float v10 {1.f};
        float v11 {0};
        float v12 {0};
        float v13 {0};
        float v14 {0};
        float v15 {1.f};
    };
    
    class Matrix4x4
    {
    public  :
        glm::mat4 m_data;

    public:
        /** Default constructor.
        @note
        It does <b>NOT</b> initialize the matrix for efficiency.
        */
        Matrix4x4() : m_data(1.0f) {} // 单位矩阵
        
        Matrix4x4(const Matrix4x4_& mat)
        {
            m_data[0][0] = mat.v0;
            m_data[1][0] = mat.v1;
            m_data[2][0] = mat.v2;
            m_data[3][0] = mat.v3;
            m_data[0][1] = mat.v4;
            m_data[1][1] = mat.v5;
            m_data[2][1] = mat.v6;
            m_data[3][1] = mat.v7;
            m_data[0][2] = mat.v8;
            m_data[1][2] = mat.v9;
            m_data[2][2] = mat.v10;
            m_data[3][2] = mat.v11;
            m_data[0][3] = mat.v12;
            m_data[1][3] = mat.v13;
            m_data[2][3] = mat.v14;
            m_data[3][3] = mat.v15;
        }

        Matrix4x4_ toMatrix4x4_() const
        {
            Matrix4x4_ res;
            res.v0  = m_data[0][0];
            res.v1  = m_data[1][0];
            res.v2  = m_data[2][0];
            res.v3  = m_data[3][0];
            res.v4  = m_data[0][1];
            res.v5  = m_data[1][1];
            res.v6  = m_data[2][1];
            res.v7  = m_data[3][1];
            res.v8  = m_data[0][2];
            res.v9  = m_data[1][2];
            res.v10 = m_data[2][2];
            res.v11 = m_data[3][2];
            res.v12 = m_data[0][3];
            res.v13 = m_data[1][3];
            res.v14 = m_data[2][3];
            res.v15 = m_data[3][3];
            return res;
        }

        Matrix4x4(const float (&float_array)[16])
        {
            m_data[0][0] = float_array[0];
            m_data[1][0] = float_array[1];
            m_data[2][0] = float_array[2];
            m_data[3][0] = float_array[3];
            m_data[0][1] = float_array[4];
            m_data[1][1] = float_array[5];
            m_data[2][1] = float_array[6];
            m_data[3][1] = float_array[7];
            m_data[0][2] = float_array[8];
            m_data[1][2] = float_array[9];
            m_data[2][2] = float_array[10];
            m_data[3][2] = float_array[11];
            m_data[0][3] = float_array[12];
            m_data[1][3] = float_array[13];
            m_data[2][3] = float_array[14];
            m_data[3][3] = float_array[15];
        }

        Matrix4x4(float m00, float m01, float m02, float m03,
                  float m10, float m11, float m12, float m13,
                  float m20, float m21, float m22, float m23,
                  float m30, float m31, float m32, float m33)
        {
            m_data[0][0] = m00; m_data[1][0] = m01; m_data[2][0] = m02; m_data[3][0] = m03;
            m_data[0][1] = m10; m_data[1][1] = m11; m_data[2][1] = m12; m_data[3][1] = m13;
            m_data[0][2] = m20; m_data[1][2] = m21; m_data[2][2] = m22; m_data[3][2] = m23;
            m_data[0][3] = m30; m_data[1][3] = m31; m_data[2][3] = m32; m_data[3][3] = m33;
        }

        Matrix4x4(const Vector4& row0, const Vector4& row1, const Vector4& row2, const Vector4& row3)
        {
            m_data[0][0] = row0.getX(); m_data[1][0] = row0.getY(); m_data[2][0] = row0.getZ(); m_data[3][0] = row0.getW();
            m_data[0][1] = row1.getX(); m_data[1][1] = row1.getY(); m_data[2][1] = row1.getZ(); m_data[3][1] = row1.getW();
            m_data[0][2] = row2.getX(); m_data[1][2] = row2.getY(); m_data[2][2] = row2.getZ(); m_data[3][2] = row2.getW();
            m_data[0][3] = row3.getX(); m_data[1][3] = row3.getY(); m_data[2][3] = row3.getZ(); m_data[3][3] = row3.getW();
        }

        Matrix4x4(const Vector3& position, const Vector3& scale, const Quaternion& rotation)
        {
            makeTransform(position, scale, rotation);
        }

        explicit Matrix4x4(const glm::mat4& mat) : m_data(mat) {}

        void fromData(const float (&float_array)[16])
        {
            m_data[0][0] = float_array[0];
            m_data[1][0] = float_array[1];
            m_data[2][0] = float_array[2];
            m_data[3][0] = float_array[3];
            m_data[0][1] = float_array[4];
            m_data[1][1] = float_array[5];
            m_data[2][1] = float_array[6];
            m_data[3][1] = float_array[7];
            m_data[0][2] = float_array[8];
            m_data[1][2] = float_array[9];
            m_data[2][2] = float_array[10];
            m_data[3][2] = float_array[11];
            m_data[0][3] = float_array[12];
            m_data[1][3] = float_array[13];
            m_data[2][3] = float_array[14];
            m_data[3][3] = float_array[15];
        }

        void toData(float (&float_array)[16]) const
        {
            float_array[0]  = m_data[0][0];
            float_array[1]  = m_data[1][0];
            float_array[2]  = m_data[2][0];
            float_array[3]  = m_data[3][0];
            float_array[4]  = m_data[0][1];
            float_array[5]  = m_data[1][1];
            float_array[6]  = m_data[2][1];
            float_array[7]  = m_data[3][1];
            float_array[8]  = m_data[0][2];
            float_array[9]  = m_data[1][2];
            float_array[10] = m_data[2][2];
            float_array[11] = m_data[3][2];
            float_array[12] = m_data[0][3];
            float_array[13] = m_data[1][3];
            float_array[14] = m_data[2][3];
            float_array[15] = m_data[3][3];
        }

        /** Creates a standard 4x4 transformation matrix with a zero translation part from a rotation/scaling 3x3
         * matrix.
         */
        void setMatrix3x3(const Matrix3x3& mat3)
        {
            // 从Matrix3x3获取数据
            m_data[0][0] = mat3[0][0]; m_data[1][0] = mat3[0][1]; m_data[2][0] = mat3[0][2]; m_data[3][0] = 0;
            m_data[0][1] = mat3[1][0]; m_data[1][1] = mat3[1][1]; m_data[2][1] = mat3[1][2]; m_data[3][1] = 0;
            m_data[0][2] = mat3[2][0]; m_data[1][2] = mat3[2][1]; m_data[2][2] = mat3[2][2]; m_data[3][2] = 0;
            m_data[0][3] = 0; m_data[1][3] = 0; m_data[2][3] = 0; m_data[3][3] = 1;
        }

        /** Creates a standard 4x4 transformation matrix with a zero translation part from a rotation/scaling
         * Quaternion.
         */
        Matrix4x4(const Quaternion& rot)
        {
            Matrix3x3 m3x3;
            rot.toRotationMatrix(m3x3);
            m_data = glm::mat4(1.0f);
            setMatrix3x3(m3x3);
        }

        // 访问元素
        float* operator[](size_t row_index)
        {
            assert(row_index < 4);
            return &m_data[0][row_index];
        }

        const float* operator[](size_t row_index) const
        {
            assert(row_index < 4);
            return &m_data[0][row_index];
        }

        // 矩阵乘法（使用GLM）
        Matrix4x4 concatenate(const Matrix4x4& m2) const
        {
            return Matrix4x4(m_data * m2.m_data);
        }

        /** Matrix concatenation using '*'.
         */
        Matrix4x4 operator*(const Matrix4x4& m2) const { return concatenate(m2); }

        /** Vector transformation using '*'.
        @remarks
        Transforms the given 3-D vector by the matrix, projecting the
        result back into <i>w</i> = 1.
        @note
        This means that the initial <i>w</i> is considered to be 1.0,
        and then all the three elements of the resulting 3-D vector are
        divided by the resulting <i>w</i>.
        */
        Vector3 operator*(const Vector3& v) const
        {
            glm::vec4 result = m_data * glm::vec4(v.getX(), v.getY(), v.getZ(), 1.0f);
            if (result.w != 0.0f) {
                result /= result.w;
            }
            return Vector3(result.x, result.y, result.z);
        }

        Vector4 operator*(const Vector4& v) const
        {
            glm::vec4 result = m_data * v.toGlm();
            return Vector4(result);
        }

        /** Matrix addition.
         */
        Matrix4x4 operator+(const Matrix4x4& m2) const
        {
            return Matrix4x4(m_data + m2.m_data);
        }

        /** Matrix subtraction.
         */
        Matrix4x4 operator-(const Matrix4x4& m2) const
        {
            return Matrix4x4(m_data - m2.m_data);
        }

        Matrix4x4 operator*(float scalar) const
        {
            return Matrix4x4(m_data * scalar);
        }

        /** Tests 2 matrices for equality.
         */
        bool operator==(const Matrix4x4& m2) const
        {
            return m_data == m2.m_data;
        }

        /** Tests 2 matrices for inequality.
         */
        bool operator!=(const Matrix4x4& m2) const
        {
            return m_data != m2.m_data;
        }

        Matrix4x4 transpose() const
        {
            return Matrix4x4(glm::transpose(m_data));
        }

        //-----------------------------------------------------------------------
        float getMinor(size_t r0, size_t r1, size_t r2, size_t c0, size_t c1, size_t c2) const
        {
            // 使用GLM的行列式计算
            glm::mat3 submatrix(
                m_data[c0][r0], m_data[c1][r0], m_data[c2][r0],
                m_data[c0][r1], m_data[c1][r1], m_data[c2][r1],
                m_data[c0][r2], m_data[c1][r2], m_data[c2][r2]
            );
            return glm::determinant(submatrix);
        }

        /*
        -----------------------------------------------------------------------
        Translation Transformation
        -----------------------------------------------------------------------
        */
        /** Sets the translation transformation part of the matrix.
         */
        void setTrans(const Vector3& v)
        {
            m_data[3][0] = v.getX();
            m_data[3][1] = v.getY();
            m_data[3][2] = v.getZ();
        }

        /** Extracts the translation transformation part of the matrix.
         */
        Vector3 getTrans() const { return Vector3(m_data[3][0], m_data[3][1], m_data[3][2]); }

        Matrix4x4 buildViewportMatrix(uint32_t width, uint32_t height)
        {
            return Matrix4x4(0.5f * (float)width,
                             0.0f,
                             0.0f,
                             0.5f * (float)width,
                             0.0f,
                             -0.5f * (float)height,
                             0.0f,
                             0.5f * (float)height,
                             0.0f,
                             0.0f,
                             -1.0f,
                             1.0f,
                             0.0f,
                             0.0f,
                             0.0f,
                             1.0f);
        }

        static Matrix4x4 mirrorMatrix(Vector4 mirror_plane)
        {
            Matrix4x4 result;
            result.m_data[0][0] = -2 * mirror_plane.getX() * mirror_plane.getX() + 1;
            result.m_data[1][0] = -2 * mirror_plane.getX() * mirror_plane.getY();
            result.m_data[2][0] = -2 * mirror_plane.getX() * mirror_plane.getZ();
            result.m_data[3][0] = 0;

            result.m_data[0][1] = -2 * mirror_plane.getY() * mirror_plane.getX();
            result.m_data[1][1] = -2 * mirror_plane.getY() * mirror_plane.getY() + 1;
            result.m_data[2][1] = -2 * mirror_plane.getY() * mirror_plane.getZ();
            result.m_data[3][1] = 0;

            result.m_data[0][2] = -2 * mirror_plane.getZ() * mirror_plane.getX();
            result.m_data[1][2] = -2 * mirror_plane.getZ() * mirror_plane.getY();
            result.m_data[2][2] = -2 * mirror_plane.getZ() * mirror_plane.getZ() + 1;
            result.m_data[3][2] = 0;

            result.m_data[0][3] = -2 * mirror_plane.getW() * mirror_plane.getX();
            result.m_data[1][3] = -2 * mirror_plane.getW() * mirror_plane.getY();
            result.m_data[2][3] = -2 * mirror_plane.getW() * mirror_plane.getZ();
            result.m_data[3][3] = 1;

            return result;
        }

        static Matrix4x4 rotationMatrix(Vector3 normal)
        {
            Vector3 up = Vector3(0, 0, 1);
            if (fabs(normal.getZ()) > 0.999f)
            {
                up = Vector3(0, 1, 0);
            }

            Vector3 left = up.cross(normal);
            up           = normal.cross(left);

            left.normalise();
            up.normalise();

            Matrix4x4 result = Matrix4x4::IDENTITY;
            result.setMatrix3x3(Matrix3x3(left, up, normal));

            return result.transpose();
        }

        /** Builds a translation matrix
         */
        void makeTrans(const Vector3& v)
        {
            m_data = glm::translate(glm::mat4(1.0f), glm::vec3(v.getX(), v.getY(), v.getZ()));
        }

        void makeTrans(float tx, float ty, float tz)
        {
            m_data = glm::translate(glm::mat4(1.0f), glm::vec3(tx, ty, tz));
        }

        /** Gets a translation matrix.
         */
        static Matrix4x4 getTrans(const Vector3& v)
        {
            return Matrix4x4(glm::translate(glm::mat4(1.0f), glm::vec3(v.getX(), v.getY(), v.getZ())));
        }

        /** Gets a translation matrix - variation for not using a vector.
         */
        static Matrix4x4 getTrans(float t_x, float t_y, float t_z)
        {
            return Matrix4x4(glm::translate(glm::mat4(1.0f), glm::vec3(t_x, t_y, t_z)));
        }

        /*
        -----------------------------------------------------------------------
        Scale Transformation
        -----------------------------------------------------------------------
        */
        /** Sets the scale part of the matrix.
         */
        void setScale(const Vector3& v)
        {
            m_data[0][0] = v.getX();
            m_data[1][1] = v.getY();
            m_data[2][2] = v.getZ();
        }

        /** Gets a scale matrix.
         */
        static Matrix4x4 getScale(const Vector3& v)
        {
            return Matrix4x4(glm::scale(glm::mat4(1.0f), glm::vec3(v.getX(), v.getY(), v.getZ())));
        }

        /** Gets a scale matrix - variation for not using a vector.
         */
        static Matrix4x4 buildScaleMatrix(float s_x, float s_y, float s_z)
        {
            return Matrix4x4(glm::scale(glm::mat4(1.0f), glm::vec3(s_x, s_y, s_z)));
        }

        /** Extracts the rotation / scaling part of the Matrix as a 3x3 matrix.
        @param m3x3 Destination Matrix3
        */
        void extract3x3Matrix(Matrix3x3& m3x3) const
        {
            m3x3[0][0] = m_data[0][0];
            m3x3[0][1] = m_data[1][0];
            m3x3[0][2] = m_data[2][0];
            m3x3[1][0] = m_data[0][1];
            m3x3[1][1] = m_data[1][1];
            m3x3[1][2] = m_data[2][1];
            m3x3[2][0] = m_data[0][2];
            m3x3[2][1] = m_data[1][2];
            m3x3[2][2] = m_data[2][2];
        }

        void extractAxes(Vector3& out_x, Vector3& out_y, Vector3& out_z) const
        {
            out_x = Vector3(m_data[0][0], m_data[0][1], m_data[0][2]);
            out_x.normalise();
            out_y = Vector3(m_data[1][0], m_data[1][1], m_data[1][2]);
            out_y.normalise();
            out_z = Vector3(m_data[2][0], m_data[2][1], m_data[2][2]);
            out_z.normalise();
        }

        /** Determines if this matrix involves a scaling. */
        bool hasScale() const
        {
            // check magnitude of column vectors (==local axes)
            float t = m_data[0][0] * m_data[0][0] + m_data[0][1] * m_data[0][1] + m_data[0][2] * m_data[0][2];
            if (!Math::realEqual(t, 1.0, (float)1e-04))
                return true;
            t = m_data[1][0] * m_data[1][0] + m_data[1][1] * m_data[1][1] + m_data[1][2] * m_data[1][2];
            if (!Math::realEqual(t, 1.0, (float)1e-04))
                return true;
            t = m_data[2][0] * m_data[2][0] + m_data[2][1] * m_data[2][1] + m_data[2][2] * m_data[2][2];
            return !Math::realEqual(t, 1.0, (float)1e-04);
        }

        /** Determines if this matrix involves a negative scaling. */
        bool hasNegativeScale() const { return determinant() < 0; }

        /** Extracts the rotation / scaling part as a quaternion from the Matrix.
         */
        Quaternion extractQuaternion() const
        {
            Matrix3x3 m3x3;
            extract3x3Matrix(m3x3);
            return Quaternion(m3x3);
        }

        Matrix4x4 adjoint() const
        {
            // 伴随矩阵 = 转置(余子式矩阵)
            float m00 = m_data[0][0], m01 = m_data[1][0], m02 = m_data[2][0], m03 = m_data[3][0];
            float m10 = m_data[0][1], m11 = m_data[1][1], m12 = m_data[2][1], m13 = m_data[3][1];
            float m20 = m_data[0][2], m21 = m_data[1][2], m22 = m_data[2][2], m23 = m_data[3][2];
            float m30 = m_data[0][3], m31 = m_data[1][3], m32 = m_data[2][3], m33 = m_data[3][3];

            float a00 =  (m11 * (m22 * m33 - m23 * m32) - m12 * (m21 * m33 - m23 * m31) + m13 * (m21 * m32 - m22 * m31));
            float a01 = -(m10 * (m22 * m33 - m23 * m32) - m12 * (m20 * m33 - m23 * m30) + m13 * (m20 * m32 - m22 * m30));
            float a02 =  (m10 * (m21 * m33 - m23 * m31) - m11 * (m20 * m33 - m23 * m30) + m13 * (m20 * m31 - m21 * m30));
            float a03 = -(m10 * (m21 * m32 - m22 * m31) - m11 * (m20 * m32 - m22 * m30) + m12 * (m20 * m31 - m21 * m30));

            float a10 = -(m01 * (m22 * m33 - m23 * m32) - m02 * (m21 * m33 - m23 * m31) + m03 * (m21 * m32 - m22 * m31));
            float a11 =  (m00 * (m22 * m33 - m23 * m32) - m02 * (m20 * m33 - m23 * m30) + m03 * (m20 * m32 - m22 * m30));
            float a12 = -(m00 * (m21 * m33 - m23 * m31) - m01 * (m20 * m33 - m23 * m30) + m03 * (m20 * m31 - m21 * m30));
            float a13 =  (m00 * (m21 * m32 - m22 * m31) - m01 * (m20 * m32 - m22 * m30) + m02 * (m20 * m31 - m21 * m30));

            float a20 =  (m01 * (m12 * m33 - m13 * m32) - m02 * (m11 * m33 - m13 * m31) + m03 * (m11 * m32 - m12 * m31));
            float a21 = -(m00 * (m12 * m33 - m13 * m32) - m02 * (m10 * m33 - m13 * m30) + m03 * (m10 * m32 - m12 * m30));
            float a22 =  (m00 * (m11 * m33 - m13 * m31) - m01 * (m10 * m33 - m13 * m30) + m03 * (m10 * m31 - m11 * m30));
            float a23 = -(m00 * (m11 * m32 - m12 * m31) - m01 * (m10 * m32 - m12 * m30) + m02 * (m10 * m31 - m11 * m30));

            float a30 = -(m01 * (m12 * m23 - m13 * m22) - m02 * (m11 * m23 - m13 * m21) + m03 * (m11 * m22 - m12 * m21));
            float a31 =  (m00 * (m12 * m23 - m13 * m22) - m02 * (m10 * m23 - m13 * m20) + m03 * (m10 * m22 - m12 * m20));
            float a32 = -(m00 * (m11 * m23 - m13 * m21) - m01 * (m10 * m23 - m13 * m20) + m03 * (m10 * m21 - m11 * m20));
            float a33 =  (m00 * (m11 * m22 - m12 * m21) - m01 * (m10 * m22 - m12 * m20) + m02 * (m10 * m21 - m11 * m20));

            return Matrix4x4(a00, a01, a02, a03,
                             a10, a11, a12, a13,
                             a20, a21, a22, a23,
                             a30, a31, a32, a33);
        }

        float determinant() const
        {
            return glm::determinant(m_data);
        }

        /** Building a Matrix4 from orientation / scale / position.
        @remarks
        Transform is performed in the order scale, rotate, translation, i.e. translation is independent
        of orientation axes, scale does not affect size of translation, rotation and scaling are always
        centered on the origin.
        */
        void makeTransform(const Vector3& position, const Vector3& scale, const Quaternion& orientation)
        {
            // 使用GLM构建变换矩阵：先缩放，再旋转，最后平移
            glm::mat4 translation = glm::translate(glm::mat4(1.0f), 
                glm::vec3(position.getX(), position.getY(), position.getZ()));
            glm::mat4 rotation = glm::mat4_cast(orientation.toGlm());
            glm::mat4 scaling = glm::scale(glm::mat4(1.0f), 
                glm::vec3(scale.getX(), scale.getY(), scale.getZ()));
            
            // 注意：GLM的矩阵乘法顺序是反的（列主序）
            // 我们需要的是 M = T * R * S
            m_data = translation * rotation * scaling;
        }

        /** Building an inverse Matrix4 from orientation / scale / position.
        @remarks
        As makeTransform except it build the inverse given the same data as makeTransform, so
        performing -translation, -rotate, 1/scale in that order.
        */
        void makeInverseTransform(const Vector3& position, const Vector3& scale, const Quaternion& orientation)
        {
            // 逆变换：先逆平移，再逆旋转，最后逆缩放
            Vector3 invScale(1.0f / scale.getX(), 1.0f / scale.getY(), 1.0f / scale.getZ());
            Quaternion invRotation = orientation.inverse();
            Vector3 invPosition = -position;
            
            makeTransform(invPosition, invScale, invRotation);
        }

        /** Decompose a Matrix4 to orientation / scale / position.
         */
        void decomposition(Vector3& position, Vector3& scale, Quaternion& orientation) const
        {
            // 使用GLM的矩阵分解
            glm::vec3 glmScale, glmTranslation, glmSkew;
            glm::quat glmRotation;
            glm::vec4 glmPerspective;
            
            glm::decompose(m_data, glmScale, glmRotation, glmTranslation, glmSkew, glmPerspective);
            
            position = Vector3(glmTranslation.x, glmTranslation.y, glmTranslation.z);
            scale = Vector3(glmScale.x, glmScale.y, glmScale.z);
            orientation = Quaternion(glmRotation);
        }

        void decompositionWithoutScale(Vector3& position, Quaternion& rotation) const
        {
            Vector3 scale;
            decomposition(position, scale, rotation);
        }

        /** Check whether or not the matrix is affine matrix.
        @remarks
        An affine matrix is a 4x4 matrix with row 3 equal to (0, 0, 0, 1),
        e.g. no projective coefficients.
        */
        bool isAffine(void) const
        {
            return m_data[0][3] == 0 && m_data[1][3] == 0 && m_data[2][3] == 0 && m_data[3][3] == 1;
        }

        /** Returns the inverse of the affine matrix.
        @note
        The matrix must be an affine matrix. @see Matrix4::isAffine.
        */
        Matrix4x4 inverseAffine() const
        {
            assert(isAffine());
            // GLM没有直接的affineInverse，我们手动计算
            // 对于仿射矩阵 M = [R t; 0 1]，逆矩阵为 [R^T -R^T*t; 0 1]
            Matrix4x4 result;
            
            // 提取旋转部分
            glm::mat3 rotation = glm::mat3(m_data);
            glm::mat3 rotationTranspose = glm::transpose(rotation);
            
            // 提取平移部分
            glm::vec3 translation = glm::vec3(m_data[3][0], m_data[3][1], m_data[3][2]);
            
            // 计算逆平移
            glm::vec3 invTranslation = -rotationTranspose * translation;
            
            // 构建逆矩阵
            result.m_data = glm::mat4(rotationTranspose);
            result.m_data[3][0] = invTranslation.x;
            result.m_data[3][1] = invTranslation.y;
            result.m_data[3][2] = invTranslation.z;
            result.m_data[0][3] = 0;
            result.m_data[1][3] = 0;
            result.m_data[2][3] = 0;
            result.m_data[3][3] = 1;
            
            return result;
        }

        /** Concatenate two affine matrices.
        @note
        The matrices must be affine matrix. @see Matrix4::isAffine.
        */
        Matrix4x4 concatenateAffine(const Matrix4x4& m2) const
        {
            assert(isAffine() && m2.isAffine());
            return Matrix4x4(m_data * m2.m_data);
        }

        /** 3-D Vector transformation specially for an affine matrix.
        @remarks
        Transforms the given 3-D vector by the matrix, projecting the
        result back into <i>w</i> = 1.
        @note
        The matrix must be an affine matrix. @see Matrix4::isAffine.
        */
        Vector3 transformAffine(const Vector3& v) const
        {
            assert(isAffine());
            glm::vec4 result = m_data * glm::vec4(v.getX(), v.getY(), v.getZ(), 1.0f);
            return Vector3(result.x, result.y, result.z);
        }

        /** 4-D Vector transformation specially for an affine matrix.
        @note
        The matrix must be an affine matrix. @see Matrix4::isAffine.
        */
        Vector4 transformAffine(const Vector4& v) const
        {
            assert(isAffine());
            glm::vec4 result = m_data * v.toGlm();
            return Vector4(result);
        }

        Matrix4x4 inverse() const
        {
            return Matrix4x4(glm::inverse(m_data));
        }

        Vector3 transformCoord(const Vector3& v)
        {
            Vector4 temp(v, 1.0f);
            Vector4 ret = (*this) * temp;
            if (ret.getW() == 0.0f)
            {
                return Vector3::ZERO;
            }
            else
            {
                ret = ret / ret.getW();
                return Vector3(ret.getX(), ret.getY(), ret.getZ());
            }
        }

        // 转换为GLM类型
        const glm::mat4& toGlm() const { return m_data; }
        glm::mat4& toGlm() { return m_data; }

        // 静态常量
        static const Matrix4x4 ZERO;
        static const Matrix4x4 ZEROAFFINE;
        static const Matrix4x4 IDENTITY;
    };

    // 全局运算符
    inline Vector4 operator*(const Vector4& v, const Matrix4x4& mat)
    {
        return mat * v;
    }

    // 静态常量定义
    inline const Matrix4x4 Matrix4x4::ZERO = Matrix4x4(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    inline const Matrix4x4 Matrix4x4::ZEROAFFINE = Matrix4x4(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1);
    inline const Matrix4x4 Matrix4x4::IDENTITY = Matrix4x4(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1);
} // namespace Lizeral

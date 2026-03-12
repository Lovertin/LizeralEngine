#pragma once
#include "runtime/function/ecs/components/component.h"
#include "runtime/core/math/matrix4.h"
#include "runtime/core/math/vector2.h"
#include "runtime/function/physics/phsicsEntityheaders.h"
#include <math.h>

// define dirty digits
const uint32_t CAMERA_DIRTY_FOV   = 1 << 1; // fov change
const uint32_t CAMERA_DIRTY_ASPECT   = 1 << 2; // apsect change
const uint32_t CAMERA_DIRTY_ZNEAR = 1 << 3; // zNear change
const uint32_t CAMERA_DIRTY_ZFAR = 1 << 4; //zFar change
const uint32_t CAMERA_DIRTY_MAIN = 1 << 5; //maincamera change

namespace Lizeral{
    REFLECTION_TYPE(CameraComponent)
    CLASS(CameraComponent:public Component,WhiteListFields){
        REFLECTION_BODY(CameraComponent)

    public:
        BEGIN_REFLECTION_UPDATED()
            REFLECTION_BIND_DIRTY("m_fov",CAMERA_DIRTY_FOV)
            REFLECTION_BIND_DIRTY("m_aspect",CAMERA_DIRTY_ASPECT)
            REFLECTION_BIND_DIRTY("m_zNear",CAMERA_DIRTY_ZNEAR)
            REFLECTION_BIND_DIRTY("m_zFar",CAMERA_DIRTY_ZFAR)
        END_REFLECTION_UPDATED()
    
        //setter
        void setFov(float fov){
            if(m_fov!=fov){ m_fov=fov ; setDirty(CAMERA_DIRTY_FOV) ;}
        }

        void setAspect(float aspect){
            if(m_aspect!=aspect){ m_aspect=aspect; setDirty(CAMERA_DIRTY_ASPECT); }
        }

        void setzNear(float zNear){
            if(m_zNear!=zNear) {m_zNear=zNear; setDirty(CAMERA_DIRTY_ZNEAR); }
        }

        void setzFar(float zFar){
            if(m_zFar!=zFar){ m_zFar=zFar; setDirty(CAMERA_DIRTY_ZFAR); }
        }

        void setMain(bool main){
            if(isMainCamera!=main) { isMainCamera=main; setDirty(CAMERA_DIRTY_MAIN); }
        }

        //without reflection
        void setViewMatrix(const Matrix4x4& mat) { m_viewMatrix = mat; }

        void setProjectionMatrix(const Matrix4x4& mat) { m_projMatrix = mat; }


        //getter
        float getFov()const {return m_fov; }

        float getAspect()const { return m_aspect; }

        float getzNear()const { return m_zNear; }

        float getzFar()const { return m_zFar; }

        const Matrix4x4& getViewMatrix() const { return m_viewMatrix; }

        const Matrix4x4& getProjectionMatrix() const { return m_projMatrix; }

        bool isMain() const {return isMainCamera;}

        //build perspective
        Matrix4x4 buildPerspective(float fov, float aspect, float zNear, float zFar) {
            float tanHalfFov = tan(fov * 0.5f);
            
            Matrix4x4 result = Matrix4x4::ZERO;
            result.m_mat[0][0] = 1.0f / (aspect * tanHalfFov);
            result.m_mat[1][1] = 1.0f / tanHalfFov;
            
            // 正确的 Z 轴深度映射
            result.m_mat[2][2] = -(zFar + zNear) / (zFar - zNear);
            result.m_mat[2][3] = -(2.0f * zFar * zNear) / (zFar - zNear);
            
            // 极其关键：W = -Z。因为我们看向 -Z 轴，负负得正，保证前方的物体 W > 0
            result.m_mat[3][2] = -1.0f; 
            result.m_mat[3][3] = 0.0f;
            
            return result;
        }

        Ray ScreenPointToRay(const Vector2& screenPos, const Vector2& screenSize) const {
            // 1. 归一化设备坐标 (NDC) [-1, 1]
            // 注意：Y轴通常需要翻转，因为屏幕原点在左上，OpenGL NDC原点在中心且Y向上
            float x = (2.0f * screenPos.x) / screenSize.x - 1.0f;
            float y = 1.0f - (2.0f * screenPos.y) / screenSize.y; 

            // 2. 裁剪空间坐标 (Clip Space)
            // 射线的起点在近平面 (z=-1)，终点在远平面 (z=1)
            // 但其实我们只需要方向。这里有一种更简单的做法：
            // 在 View Space 中，相机原点是 (0,0,0)。
            // 我们只需要算出鼠标点击点在 View Space 中的方向。

            // 逆投影矩阵
            Matrix4x4 invProj = m_projMatrix.inverse();
            Matrix4x4 invView = m_viewMatrix.inverse(); // CameraComponent 应该缓存了这个

            Vector4 rayClip(x, y, -1.0f, 1.0f);
            Vector4 rayEye = invProj * rayClip;
            
            // 我们只需要方向，所以设 z=-1 (指向前方), w=0
            rayEye = Vector4(rayEye.x, rayEye.y, -1.0f, 0.0f);

            // 3. 转换到世界空间
            Vector3 rayDir = (invView * rayEye).toVector3();
            rayDir.normalise();

            // 4. 射线的起点就是相机的世界坐标
            // 注意：如果是正交相机(2D)，起点不一样。这里假设是透视相机。
            // CameraComponent 本身不存位置，位置在 Transform 里，这里可能需要传入 Transform 或者 ViewMatrix 的逆矩阵的平移部分
            Vector3 origin = Vector3(invView[0][3], invView[1][3], invView[2][3]); 

            return Ray(origin, rayDir);
        }

        // 构建无限远反向深度的透视投影矩阵 (假设为右手坐标系，相机看向 -Z)
        Matrix4x4 BuildPerspectiveInfiniteReverseZ(float fov, float aspect, float zNear) {
            Matrix4x4 result; 
            float tanHalfFovy = std::tan(fov * 0.5f);
            
            result[0][0] = 1.0f / (aspect * tanHalfFovy);
            result[1][1] = 1.0f / tanHalfFovy; 

            result[2][2] = 0.0f;  
            result[2][3] = zNear; 

            result[3][2] = -1.0f; 
            result[3][3] = 0.0f;  
            
            return result;
        }
    
    private:
        META(Enable)
        float m_fov {60.0f};

        META(Enable)
        float m_aspect {0.75f};

        META(Enable)
        float m_zNear {0.1f};

        META(Enable)
        float m_zFar {1000.0f};

        Matrix4x4 m_viewMatrix;
        Matrix4x4 m_projMatrix;

        META(Enable)
        bool isMainCamera {true};
    
    };
}
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

            result.m_mat[2][2] = -(zFar + zNear) / (zFar - zNear);
            result.m_mat[2][3] = -(2.0f * zFar * zNear) / (zFar - zNear);

            result.m_mat[3][2] = -1.0f; 
            result.m_mat[3][3] = 0.0f;
            
            return result;
        }

        Ray ScreenPointToRay(const Vector2& screenPos, const Vector2& screenSize) const {

            float x = (2.0f * screenPos.x) / screenSize.x - 1.0f;
            float y = 1.0f - (2.0f * screenPos.y) / screenSize.y; 

            Matrix4x4 invProj = m_projMatrix.inverse();
            Matrix4x4 invView = m_viewMatrix.inverse();

            Vector4 rayClip(x, y, -1.0f, 1.0f);
            Vector4 rayEye = invProj * rayClip;

            rayEye = Vector4(rayEye.x, rayEye.y, -1.0f, 0.0f);

            Vector3 rayDir = (invView * rayEye).toVector3();
            rayDir.normalise();

            Vector3 origin = Vector3(invView[0][3], invView[1][3], invView[2][3]); 

            return Ray(origin, rayDir);
        }

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
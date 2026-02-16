#pragma once
#include "runtime/function/ecs/components/component.h"
#include "runtime/core/math/matrix4.h"

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
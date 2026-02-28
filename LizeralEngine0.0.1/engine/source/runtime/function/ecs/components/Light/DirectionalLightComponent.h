#pragma once
#include "runtime/core/meta/reflection/reflection.h"
#include "runtime/function/ecs/components/component.h"
#include "runtime/core/math/vector3.h"

namespace Lizeral{
    REFLECTION_TYPE(DirectionLightComponent)
    CLASS(DirectionLightComponent : public Component,WhiteListFields){
        REFLECTION_BODY(DirectionLightComponent)
        public:
            //setter
            void setColor(const Vector3& color){ m_color=color; }
            void setColor(const float r,const float g,const float b){ m_color=Vector3(r,g,b); }
            void setIntensity(const float intensity){ m_intensity=intensity; }
            // void setDirection(const Vector3& direction){m_direction=direction;}
            // void setDirection(const float x,const float y,const float z){m_direction=Vector3(x,y,z);}
            void setGlobal(const bool main){ m_isGlobal=main; }

            //getter
            Vector3 getColor() const { return m_color;}
            float getIntensity() const { return m_intensity; }
            bool isGlobal() const { return m_isGlobal; }


  
        private:
            META(Enable)
            Vector3 m_color{1.0f,1.0f,1.0f};

            META(Enable)
            float m_intensity{1.0f};
            // Vector3 m_direction{0.5f, -1.0f, 0.5f};

            META(Enable)
            bool m_isGlobal{true};
    };

}
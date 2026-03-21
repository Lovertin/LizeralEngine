#pragma once
#include "runtime/function/ecs/components/component.h"
#include "runtime/core/math/vector3.h"

namespace Lizeral{

    REFLECTION_TYPE(PointLightComponent)
    CLASS(PointLightComponent:public Component,WhiteListFields){
        REFLECTION_BODY(PointLightComponent)
    public:
        void setSpotIntensity(const float intensity) { m_spotintensity = intensity ;}
        void setSpotLightColor(const Vector3& color) { m_spotlightcolor = color ;}
        void setSpotRadius(const float radius) { m_spotradius = radius ;}

        float getSpotIntensity() const { return m_spotintensity ;}
        Vector3 getSpotLightColor() const { return m_spotlightcolor ;}
        float getSpotRadius() const { return m_spotradius ;}

    private:
        META(Enable)
        float m_spotintensity{1.0f};

        META(Enable)
        Vector3 m_spotlightcolor{1.0f,1.0f,1.0f};

        META(Enable)
        float m_spotradius {5.0f};
    };
}

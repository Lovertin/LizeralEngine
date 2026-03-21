#pragma once
#include "runtime/core/math/vector3.h"

namespace Lizeral{
    struct Ray {
        Vector3 origin; 
        Vector3 direction; 

        Ray() = default;
        Ray(const Vector3& o, const Vector3& d) : origin(o), direction(d) {}

        Vector3 GetPoint(float distance) const {
            return origin + direction * distance;
        }
    };
}
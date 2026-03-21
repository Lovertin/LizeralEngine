// runtime/function/physics/raycast_hit.h
#pragma once
#include "runtime/core/math/vector3.h"
#include "runtime/function/ecs/entity.h" 

namespace Lizeral {
    struct RaycastHit {
        Entity entity;  
        Vector3 point; 
        Vector3 normal;
        float distance; 
        bool hasHit;      

        RaycastHit() : hasHit(false), distance(0.0f) {}
    };
}
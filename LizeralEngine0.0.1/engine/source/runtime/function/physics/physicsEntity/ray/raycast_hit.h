// runtime/function/physics/raycast_hit.h
#pragma once
#include "runtime/core/math/vector3.h"
#include "runtime/function/ecs/entity.h" // 假设你有 Entity ID 定义

namespace Lizeral {
    struct RaycastHit {
        Entity entity;      // 击中了哪个实体 (关键！)
        Vector3 point;      // 击中点的世界坐标
        Vector3 normal;     // 击中表面的法线 (用于做弹孔、反弹效果)
        float distance;     // 距离原点的长度
        bool hasHit;        // 是否击中

        RaycastHit() : hasHit(false), distance(0.0f) {}
    };
}
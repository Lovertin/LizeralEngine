#pragma once
#include "runtime/core/math/vector3.h"

namespace Lizeral{
    struct Ray {
        Vector3 origin;    // 起点
        Vector3 direction; // 方向 (必须归一化)

        Ray() = default;
        Ray(const Vector3& o, const Vector3& d) : origin(o), direction(d) {}

        // 辅助函数：获取射线上某一点的位置
        Vector3 GetPoint(float distance) const {
            return origin + direction * distance;
        }
    };
}
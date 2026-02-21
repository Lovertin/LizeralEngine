#pragma once

#include "runtime/function/ecs/components/component.h"
#include "runtime/function/pipeline/mesh/mesh.h"
#include "runtime/function/pipeline/shader/shader.h"

namespace Lizeral{
   struct MeshRendererComponent : public Component {
        Mesh* mesh;       // 指向模型数据
        Shader* shader;   // 指向着色器
        Vector3 color;    // 材质参数（最简单的形式）
    };
}
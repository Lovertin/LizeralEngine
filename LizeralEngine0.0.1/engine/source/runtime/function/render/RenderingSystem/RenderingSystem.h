#pragma once
#include <glad/glad.h>
#include "runtime/function/ecs/registry.h"
#include "runtime/function/res_type/shader/shader.h"
#include <iostream>

namespace Lizeral{
    class CameraComponent;
}

namespace Lizeral {

    class RenderingSystem {
    private:
        // --- 阴影管线核心资源 ---
        unsigned int m_depthMapFBO = 0; // 深度缓冲对象
        unsigned int m_depthMapArray = 0;    // 真正的阴影深度贴图
        
        const unsigned int SHADOW_WIDTH = 2048;  // 阴影贴图分辨率
        const unsigned int SHADOW_HEIGHT = 2048;

        std::shared_ptr<Shader> m_shadowShader;

        const std::vector<float> shadowCascadeLevels{ 15.0f, 50.0f, 150.0f }; // 视锥体的 3 个分割点（切蛋糕的位置）
        // const int NUM_CASCADES = shadowCascadeLevels.size(); // 3 级级联
        const int NUM_CASCADES = 4;

        unsigned int m_DefaultFBO = 0; // 默认为 0，兼容 Sandbox

    public:
        RenderingSystem() = default;
        ~RenderingSystem() { Shutdown(); }

        // 初始化渲染系统（分配 FBO 等硬件资源）
        void Initialize();

        // 核心渲染循环
        void Tick(Registry& registry, float deltaTime);

        // 清理显存
        void Shutdown();

        Matrix4x4  getLightSpaceMatrix(
            const float nearPlane, 
            const float farPlane, 
            const CameraComponent& cam, 
            const Matrix4x4& camView, 
            const Vector3& lightDir) ;

        void SetDefaultFBO(unsigned int fbo) { m_DefaultFBO = fbo; }
        
    };

} // namespace Lizeral
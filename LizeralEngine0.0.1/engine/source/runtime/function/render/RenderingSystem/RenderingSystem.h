#pragma once
#include <glad/glad.h>
#include "runtime/function/ecs/registry.h"
#include "runtime/function/res_type/shader/shader.h"
#include "runtime/function/res_type/Texture/Texture2D.h"
#include "runtime/function/physics/physicsDebug/PhysicsDebugDrawer.h"
#include "runtime/function/physics/PhysicsSystem.h"

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

        int m_viewX = 0, m_viewY = 0, m_viewW = 1920, m_viewH = 1080;

        std::shared_ptr<Shader> m_skyboxShader;
        std::shared_ptr<Texture2D> m_skyboxTexture; // 直接使用资源系统！

        unsigned int m_cubeVAO = 0;
        unsigned int m_cubeVBO = 0;

        void SetupSkyboxCube();

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

        void SetViewport(int x, int y, int w, int h) {
            m_viewX = x; m_viewY = y; m_viewW = w; m_viewH = h;
        }
        
    };

} // namespace Lizeral
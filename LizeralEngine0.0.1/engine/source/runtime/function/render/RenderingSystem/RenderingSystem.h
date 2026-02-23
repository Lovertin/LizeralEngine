#pragma once
#include <glad/glad.h>
#include "runtime/function/ecs/registry.h"
#include "runtime/function/res_type/shader/shader.h"
#include <iostream>

namespace Lizeral {

    class RenderingSystem {
    private:
        // --- 阴影管线核心资源 ---
        unsigned int m_depthMapFBO = 0; // 深度缓冲对象
        unsigned int m_depthMap = 0;    // 真正的阴影深度贴图
        
        const unsigned int SHADOW_WIDTH = 2048;  // 阴影贴图分辨率
        const unsigned int SHADOW_HEIGHT = 2048;

        std::shared_ptr<Shader> m_shadowShader;

    public:
        RenderingSystem() = default;
        ~RenderingSystem() { Shutdown(); }

        // 初始化渲染系统（分配 FBO 等硬件资源）
        void Initialize();

        // 核心渲染循环
        void Tick(Registry& registry, float deltaTime);

        // 清理显存
        void Shutdown();
    };

} // namespace Lizeral
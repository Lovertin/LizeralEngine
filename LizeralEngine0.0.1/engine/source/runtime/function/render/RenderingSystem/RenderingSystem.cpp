#include "RenderingSystem.h"
#include "runtime/function/ecs/components/componentAll.h"
#include "runtime/function/res_type/Material/PBRMaterial.h"
#include "runtime/core/math/matrix4.h"
#include "runtime/core/math/vector3.h"

namespace Lizeral {

    void RenderingSystem::Initialize() {

        m_shadowShader = std::make_shared<Shader>(
            "C:\\Lizeral Engine\\LizeralEngine0.0.1\\engine\\source\\runtime\\function\\sandbox\\shaderTest\\shadow.vert", 
            "C:\\Lizeral Engine\\LizeralEngine0.0.1\\engine\\source\\runtime\\function\\sandbox\\shaderTest\\shadow.frag"
        );
        // ==========================================
        // 【核心】：创建 Shadow Map FBO 与 深度贴图
        // ==========================================

        // 1. 生成帧缓冲对象 (FBO)
        glGenFramebuffers(1, &m_depthMapFBO);

        // 2. 生成一张纹理，用来存储深度信息
        glGenTextures(1, &m_depthMap);
        glBindTexture(GL_TEXTURE_2D, m_depthMap);
        
        // 申请一块 2048x2048 的显存，格式为 GL_DEPTH_COMPONENT (只存深度，不存颜色)
        glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT32F, SHADOW_WIDTH, SHADOW_HEIGHT, 0, GL_DEPTH_COMPONENT, GL_FLOAT, NULL);
        
       glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
        float borderColor[] = { 1.0f, 1.0f, 1.0f, 1.0f };
        glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, borderColor);

        // 3. 将这张深度纹理绑定到我们的 FBO 上
        glBindFramebuffer(GL_FRAMEBUFFER, m_depthMapFBO);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, m_depthMap, 0);
        
        // unsigned int rbo;
        // glGenRenderbuffers(1, &rbo);
        // glBindRenderbuffer(GL_RENDERBUFFER, rbo);
        // glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT32F, SHADOW_WIDTH, SHADOW_HEIGHT);
        // glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, rbo);

        // 【极其重要】：告诉 OpenGL，这个 FBO 我们不画颜色，只画深度！
        // 如果不加这两句，OpenGL 会因为找不到颜色附件而报错。
        glDrawBuffer(GL_NONE);
        glReadBuffer(GL_NONE);

        // 检查 FBO 是否创建成功
        if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
            std::cerr << "[RenderingSystem] Error: Shadow Map Framebuffer is not complete!" << std::endl;
        }

        // 解绑 FBO，恢复到默认的屏幕渲染
        glBindFramebuffer(GL_FRAMEBUFFER, 0);

        std::cout << "[RenderingSystem] Initialized. Shadow Map (2048x2048) Created." << std::endl;
    }

    void RenderingSystem::Tick(Registry& registry, float deltaTime) {
        Matrix4x4 viewMatrix = Matrix4x4::IDENTITY;
        Matrix4x4 projMatrix = Matrix4x4::IDENTITY;
        Vector3 camPos(0, 0, 0);

        auto cameraView = registry.view<CameraComponent, TransformComponent>();
        for (auto entity : cameraView) {
            auto& cam = cameraView.get<CameraComponent>(entity);
            if (cam.isMain()) {
                auto& t = cameraView.get<TransformComponent>(entity);
                viewMatrix = cam.getViewMatrix();
                projMatrix = cam.getProjectionMatrix();
                camPos = t.getPosition();
                break; // 找到主相机就退出循环
            }
        }

        Vector3 lightDir(0.0f, 1.0f, 0.001f);
        lightDir.normalise();
        Vector3 lightColor(3.0f, 3.0f, 3.0f);

        Vector3 targetPos(20.0f, 0.0f, 0.0f);

        Vector3 lightPos(100.0f, 50.0f, 40.0f);

        Matrix4x4 lightProjection = Matrix4x4::ortho(-15.0f, 15.0f, -15.0f, 15.0f, 1.0f, 150.0f);
        // 让太阳看向场景原点 (0,0,0)
        Matrix4x4 lightView = Matrix4x4::lookAt(lightPos, targetPos, Vector3(0.0f, 0.0f, 1.0f));
        
        // 光源空间的终极矩阵 (传给 Shader 的 lightSpaceMatrix)
        Matrix4x4 lightSpaceMatrix = lightProjection * lightView;

        // ==========================================
        // Pass 1: 深度渲染阶段 (Shadow Pass)
        // 目标：从太阳的视角画一遍场景，不画颜色，只把深度写入 FBO
        // ==========================================
        // 1. 切换到我们的 FBO
        glBindFramebuffer(GL_FRAMEBUFFER, m_depthMapFBO);
        // 【关键】：视口必须和阴影贴图的分辨率一模一样！
        glViewport(0, 0, SHADOW_WIDTH, SHADOW_HEIGHT); 

        glEnable(GL_DEPTH_TEST); 
        glDepthMask(GL_TRUE);

        glClear(GL_DEPTH_BUFFER_BIT); // 只需要清空深度！
        
        // 开启剔除正面，解决“阴影悬浮(Peter Panning)”的经典伪影问题
        glEnable(GL_CULL_FACE);
        glCullFace(GL_FRONT); 

        // 假设你有 m_shadowShader 这个成员变量
        m_shadowShader->Bind();
        m_shadowShader->SetUniformMat4f("lightSpaceMatrix", lightSpaceMatrix);

        auto modelView = registry.view<TransformComponent, ModelComponent>();
        for (auto entity : modelView) {
            auto& t = modelView.get<TransformComponent>(entity);
            auto& m = modelView.get<ModelComponent>(entity);
            if (!m.m_Model) continue;

            m_shadowShader->SetUniformMat4f("model", t.getMatrix());
            
            // 用极简的 shadowShader 把这个实体的所有网格画一遍
            for (auto& mesh : m.m_Model->GetMeshes()) {
                mesh.Draw(); 
            }
        }

        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        glViewport(0, 0, 1280, 720);

        glEnable(GL_DEPTH_TEST);
        glDisable(GL_CULL_FACE); // 根据你的模型情况，可以开启背面剔除
        
        for (auto entity : modelView) {
            auto& t = modelView.get<TransformComponent>(entity);
            auto& m = modelView.get<ModelComponent>(entity);

            if (!m.m_Model) continue;

            Matrix4x4 modelMatrix = t.getMatrix();

            // 遍历模型的所有子网格
            for (auto& mesh : m.m_Model->GetMeshes()) {
                
                // 将材质向下转型为 PBRMaterial (未来可以通过 Material 虚函数多态处理，无需转型)
                auto pbrMat = std::dynamic_pointer_cast<PBRMaterial>(mesh.m_Material);
                if (!pbrMat || !pbrMat->GetShader()) continue; // 防御性编程

                auto shader = pbrMat->GetShader();
                shader->Bind();

                // A. 注入引擎级/系统级数据 (所有物体都一样的环境数据)
                shader->SetUniformMat4f("model", modelMatrix);
                shader->SetUniformMat4f("view", viewMatrix);
                shader->SetUniformMat4f("projection", projMatrix);
                shader->SetUniformVector3f("u_camPos", camPos);
                
                shader->SetUniformVector3f("u_lightDir", lightDir);
                shader->SetUniformVector3f("u_lightColor", lightColor);

                shader->SetUniformMat4f("lightSpaceMatrix", lightSpaceMatrix);
                glActiveTexture(GL_TEXTURE7); // 用一个空闲的纹理槽位
                glBindTexture(GL_TEXTURE_2D, m_depthMap);
                shader->SetUniform1i("shadowMap", 7);

                // B. 注入材质级数据 (数据驱动的核心！让材质把自己的参数喂给 Shader)
                pbrMat->BindAndApply();

                // C. 提交 Draw Call
                // 注意：这里不要调用 m.m_Model->Draw()，因为那是画整个模型。我们现在是按 Mesh 拆分画的。
                // 假设你的 Mesh 类有类似 Draw() 或 BindVAO() + glDrawElements() 的方法
                mesh.Draw(); 
            }
        }

        glUseProgram(0);

        // glActiveTexture(GL_TEXTURE0);
        // // 【注意】：请确保你能拿到 m_depthMap 的 ID！
        // glBindTexture(GL_TEXTURE_2D, m_depthMap); 

        // // 切换到 2D 正交投影模式 (假设你的屏幕是 1280x720)
        // glMatrixMode(GL_PROJECTION);
        // glPushMatrix();
        // glLoadIdentity();
        // glOrtho(0, 1280, 0, 720, -1, 1);
        // glMatrixMode(GL_MODELVIEW);
        // glPushMatrix();
        // glLoadIdentity();

        // // 关闭光照和深度测试，纯粹画一张 2D 图片
        // glDisable(GL_DEPTH_TEST);
        // glDisable(GL_LIGHTING);
        // glEnable(GL_TEXTURE_2D);
        // glColor3f(1.0f, 1.0f, 1.0f); 

        // // 在右上角画一个 300x300 的正方形来显示深度图
        // glBegin(GL_QUADS);
        // glTexCoord2f(0, 0); glVertex2f(980, 420);
        // glTexCoord2f(1, 0); glVertex2f(1280, 420);
        // glTexCoord2f(1, 1); glVertex2f(1280, 720);
        // glTexCoord2f(0, 1); glVertex2f(980, 720);
        // glEnd();

        // // 恢复状态
        // glEnable(GL_DEPTH_TEST);
        // glPopMatrix();
        // glMatrixMode(GL_PROJECTION);
        // glPopMatrix();
    }

    void RenderingSystem::Shutdown() {
        if (m_depthMapFBO != 0) {
            glDeleteFramebuffers(1, &m_depthMapFBO);
            m_depthMapFBO = 0;
        }
        if (m_depthMap != 0) {
            glDeleteTextures(1, &m_depthMap);
            m_depthMap = 0;
        }
    }

} // namespace Lizeral
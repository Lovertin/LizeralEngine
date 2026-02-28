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
        // 【核心】：创建 CSM Shadow Map Array FBO
        // ==========================================

        // 1. 生成帧缓冲对象 (FBO)
        glGenFramebuffers(1, &m_depthMapFBO);

        // 2. 生成纹理数组，用来存储多层深度信息
        glGenTextures(1, &m_depthMapArray);
        // 【修复 1】：必须绑定到 GL_TEXTURE_2D_ARRAY
        glBindTexture(GL_TEXTURE_2D_ARRAY, m_depthMapArray); 
        
        // 申请显存 (长, 宽, 层数)
        glTexImage3D(GL_TEXTURE_2D_ARRAY, 0, GL_DEPTH_COMPONENT32F, 
                     SHADOW_WIDTH, SHADOW_HEIGHT, NUM_CASCADES, 
                     0, GL_DEPTH_COMPONENT, GL_FLOAT, NULL);
        
        // 【修复 2】：参数设置也必须全部对着 GL_TEXTURE_2D_ARRAY
        glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
        glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
        float borderColor[] = { 1.0f, 1.0f, 1.0f, 1.0f };
        glTexParameterfv(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_BORDER_COLOR, borderColor);

        // 3. 将深度纹理绑定到 FBO
        glBindFramebuffer(GL_FRAMEBUFFER, m_depthMapFBO);
        
        // 挂载整个纹理数组到 FBO (在 Tick 里面会用 glFramebufferTextureLayer 逐层绘制)
        glFramebufferTexture(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, m_depthMapArray, 0);

        // 告诉 OpenGL 这个 FBO 只画深度，不画颜色
        glDrawBuffer(GL_NONE);
        glReadBuffer(GL_NONE);

        // 检查 FBO 状态
        if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
            std::cerr << "[RenderingSystem] Error: Shadow Map Framebuffer is not complete!" << std::endl;
        } else {
            std::cout << "[RenderingSystem] CSM FBO Initialized Successfully. (Array Layers: " << NUM_CASCADES << ")" << std::endl;
        }

        // 解绑 FBO，恢复默认状态
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
    }

    void RenderingSystem::Tick(Registry& registry, float deltaTime) {
        // ========================================================
        // 阶段 1：搜集场景核心数据 (Data Collection)
        // ========================================================
        Matrix4x4 viewMatrix = Matrix4x4::IDENTITY;
        Matrix4x4 projMatrix = Matrix4x4::IDENTITY;
        Vector3 camPos(0, 0, 0);
        const CameraComponent* mainCamPtr = nullptr; // 缓存主相机的指针

        // 1.1 获取主摄像机
        auto cameraView = registry.view<CameraComponent, TransformComponent>();
        for (auto entity : cameraView) {
            auto& cam = cameraView.get<CameraComponent>(entity);
            if (cam.isMain()) {
                mainCamPtr = &cam;
                auto& t = cameraView.get<TransformComponent>(entity);
                viewMatrix = cam.getViewMatrix();
                projMatrix = cam.getProjectionMatrix();
                camPos = t.getPosition();
                break; 
            }
        }

        if (!mainCamPtr) return; // 如果场景里没有主相机，直接罢工，保护引擎

        // 1.2 获取全局定向光 (太阳)
        Vector3 lightDirToSun(1.0f, 1.5f, 1.0f);
        lightDirToSun.normalise();
        Vector3 lightColor(3.0f, 3.0f, 3.0f);  // 默认值

        auto sunlightView = registry.view<TransformComponent, DirectionLightComponent>();
        for (auto entity : sunlightView) {
            auto& light = sunlightView.get<DirectionLightComponent>(entity);
            if (light.isGlobal()) {
                auto& sunTransform = registry.get<TransformComponent>(entity);
                lightDirToSun = sunTransform.getForward(); // 动态获取朝向

                lightColor = light.getColor() * light.getIntensity(); // 动态获取光强
                // lightColor = Vector3(3.0f, 3.0f, 3.0f);
                // lightDirToSun = Vector3(1.5f, 1.0f, 1.5f);
                // lightDirToSun.normalise();
                break; // 只取第一个全局光
            }
        }

        // ========================================================
        // 阶段 2：CSM 数学核爆 (Calculate Cascade Matrices)
        // ========================================================
        // 这里的距离必须和你在 shadow.frag 里接收的级联距离严格对应
        std::vector<float> shadowCascadeLevels{ 15.0f, 50.0f, 150.0f }; 
        std::vector<Matrix4x4> lightSpaceMatrices;
        lightSpaceMatrices.reserve(4);

        // 利用刚刚缓存的 mainCamPtr 计算 4 层级联的包围盒矩阵
        // 注意最后一层的 farPlane 应该尽量贴合摄像机的实际 zFar，这里写 1000.0f 或 mainCamPtr->getzFar()
        lightSpaceMatrices.push_back(getLightSpaceMatrix(mainCamPtr->getzNear(), shadowCascadeLevels[0], *mainCamPtr, viewMatrix, lightDirToSun));
        lightSpaceMatrices.push_back(getLightSpaceMatrix(shadowCascadeLevels[0], shadowCascadeLevels[1], *mainCamPtr, viewMatrix, lightDirToSun));
        lightSpaceMatrices.push_back(getLightSpaceMatrix(shadowCascadeLevels[1], shadowCascadeLevels[2], *mainCamPtr, viewMatrix, lightDirToSun));
        lightSpaceMatrices.push_back(getLightSpaceMatrix(shadowCascadeLevels[2], mainCamPtr->getzFar(), *mainCamPtr, viewMatrix, lightDirToSun));


        // ========================================================
        // 阶段 3：深度渲染管线 (Shadow Pass - 连拍 4 张)
        // ========================================================
        glBindFramebuffer(GL_FRAMEBUFFER, m_depthMapFBO);
        glViewport(0, 0, SHADOW_WIDTH, SHADOW_HEIGHT); 

        glEnable(GL_DEPTH_TEST); 
        glDepthMask(GL_TRUE);
        glEnable(GL_CULL_FACE);
        glCullFace(GL_BACK); // 剔除正面防摩尔纹

        m_shadowShader->Bind();

        auto modelView = registry.view<TransformComponent, ModelComponent>();

        // 【CSM 核心循环】：一层一层地切蛋糕，并拍照
        for (int layer = 0; layer < lightSpaceMatrices.size(); ++layer) {
            // 【关键 API】：把 FBO 的画板切换到纹理数组的第 layer 层
            glFramebufferTextureLayer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, m_depthMapArray, 0, layer);
            glClear(GL_DEPTH_BUFFER_BIT); // 务必清空这一层的老数据！

            // 告诉 Shader 这张底片用哪个光照矩阵
            m_shadowShader->SetUniformMat4f("lightSpaceMatrix", lightSpaceMatrices[layer]);

            // 画出整个世界
            for (auto entity : modelView) {
                auto& t = modelView.get<TransformComponent>(entity);
                auto& m = modelView.get<ModelComponent>(entity);
                if (!m.m_Model) continue;

                m_shadowShader->SetUniformMat4f("model", t.getMatrix());
                for (auto& mesh : m.m_Model->GetMeshes()) {
                    mesh.Draw(); 
                }
            }
        }


        // ========================================================
        // 阶段 4：颜色渲染管线 (Color Pass)
        // ========================================================
        glBindFramebuffer(GL_FRAMEBUFFER, 0); // 切回屏幕
        glViewport(0, 0, 1280, 720); // 恢复屏幕分辨率

        glEnable(GL_DEPTH_TEST);
        glCullFace(GL_BACK); // 颜色渲染必须恢复背面剔除
        
        for (auto entity : modelView) {
            auto& t = modelView.get<TransformComponent>(entity);
            auto& m = modelView.get<ModelComponent>(entity);
            if (!m.m_Model) continue;

            Matrix4x4 modelMatrix = t.getMatrix();

            for (auto& mesh : m.m_Model->GetMeshes()) {
                
                auto pbrMat = std::dynamic_pointer_cast<PBRMaterial>(mesh.m_Material);
                if (!pbrMat || !pbrMat->GetShader()) continue; 

                auto shader = pbrMat->GetShader();
                shader->Bind();

                // --- 4.1 注入通用环境数据 ---
                shader->SetUniformMat4f("model", modelMatrix);
                shader->SetUniformMat4f("view", viewMatrix);
                shader->SetUniformMat4f("projection", projMatrix);
                shader->SetUniformVector3f("u_camPos", camPos);
                
                shader->SetUniformVector3f("u_lightDir", lightDirToSun);
                shader->SetUniformVector3f("u_lightColor", lightColor);

                // --- 4.2 注入 CSM 数据！ ---
                shader->SetUniform1i("u_cascadeCount", lightSpaceMatrices.size());
                
                for (size_t i = 0; i < lightSpaceMatrices.size(); ++i) {
                    shader->SetUniformMat4f("u_lightSpaceMatrices[" + std::to_string(i) + "]", lightSpaceMatrices[i]);
                    if (i < shadowCascadeLevels.size()) {
                        shader->SetUniform1f("u_cascadePlaneDistances[" + std::to_string(i) + "]", shadowCascadeLevels[i]);
                    }
                }

                // --- 4.3 绑定纹理数组 ---
                glActiveTexture(GL_TEXTURE7); 
                glBindTexture(GL_TEXTURE_2D_ARRAY, m_depthMapArray); // 极其关键：GL_TEXTURE_2D_ARRAY
                shader->SetUniform1i("shadowMap", 7);

                // --- 4.4 渲染网格 ---
                pbrMat->BindAndApply();
                mesh.Draw(); 
            }
        }

        glUseProgram(0);
    }

    void RenderingSystem::Shutdown() {
        if (m_depthMapFBO != 0) {
            glDeleteFramebuffers(1, &m_depthMapFBO);
            m_depthMapFBO = 0;
        }
        // if (m_depthMap != 0) {
        //     glDeleteTextures(1, &m_depthMap);
        //     m_depthMap = 0;
        // }
    }


    std::vector<Vector3> getFrustumCornersWorldSpace(const Matrix4x4& projMatrix, const Matrix4x4& viewMatrix) {
        // 1. 计算视图投影矩阵的逆矩阵
        Matrix4x4 viewProj = projMatrix * viewMatrix;
        Matrix4x4 invViewProj = viewProj.inverse();

        std::vector<Vector3> corners;
        corners.reserve(8);

        // 2. NDC 空间下的 8 个标准顶点
        // X, Y, Z 分别在 -1 和 1 之间
        for (unsigned int x = 0; x < 2; ++x) {
            for (unsigned int y = 0; y < 2; ++y) {
                for (unsigned int z = 0; z < 2; ++z) {
                    
                    // 将 0,1 映射到 -1.0, 1.0
                    float ptX = 2.0f * x - 1.0f;
                    float ptY = 2.0f * y - 1.0f;
                    float ptZ = 2.0f * z - 1.0f;

                    // 构造 NDC 坐标
                    Vector3 ptNDC(ptX, ptY, ptZ);

                    // 乘以逆矩阵，还原到世界空间
                    Vector3 ptWorld = invViewProj * ptNDC; 

                    // 因为我们是用逆矩阵还原透视投影，w 分量不一定为 1，必须除以 w 才能得到真实的 3D 坐标。
                    Vector4 ptWorld4 = invViewProj * Vector4(ptX, ptY, ptZ, 1.0f);
                    
                    corners.push_back(Vector3(ptWorld4.x / ptWorld4.w
                                            , ptWorld4.y / ptWorld4.w, 
                                              ptWorld4.z / ptWorld4.w));

                    // corners.push_back(ptWorld);
                }
            }
        }
        return corners;
    }

    Matrix4x4 RenderingSystem::getLightSpaceMatrix(
    const float nearPlane, 
    const float farPlane, 
    const CameraComponent& cam, 
    const Matrix4x4& camView, 
    const Vector3& lightDir) 
{
    // 1. USE THE EXACT CAMERA PROJECTION to prevent "backward/mirrored" frustums
    Matrix4x4 proj = const_cast<CameraComponent&>(cam).buildPerspective(cam.getFov(), cam.getAspect(), nearPlane, farPlane);
    std::vector<Vector3> corners = getFrustumCornersWorldSpace(proj, camView);

    // 2. Find the center of the sub-frustum
    Vector3 center(0.0f, 0.0f, 0.0f);
    for (const auto& v : corners) {
        center += v;
    }
    center /= corners.size();

    // 3. Calculate the Bounding Sphere Radius 
    // This perfectly encapsulates the frustum and stops shadows from "wobbling" when the camera rotates!
    float radius = 0.0f;
    for (const auto& v : corners) {
        float distance = (v - center).length();
        radius = std::max(radius, distance);
    }
    
    // Optional: Snap radius to texel increments here if you want absolutely 0 sub-pixel jitter, 
    // but the bounding sphere alone fixes 95% of the movement artifacts.

    // 4. Position the Sun Camera
    // We pull the sun back by a fixed 'zPullback' distance so it can see objects casting shadows from behind the camera.
    float zPullback = 300.0f;
    Vector3 lightPos = center + lightDir * zPullback;
    Vector3 upVec = (std::abs(lightDir.y) > 0.999f) ? Vector3(0.0f, 0.0f, 1.0f) : Vector3(0.0f, 1.0f, 0.0f);
    Matrix4x4 lightView = Matrix4x4::lookAt(lightPos, center, upVec);

    // 5. Build the Ortho Matrix tightly around the sphere
    // Because the light looks EXACTLY at 'center', the sphere perfectly maps to X/Y bounds from -radius to +radius.
    // ZERO padding is needed! Your resolution will be incredibly sharp.
    float minX = -radius;
    float maxX =  radius;
    float minY = -radius;
    float maxY =  radius;

    // 6. Perfect Z-Bounds
    // The center is 'zPullback' units away from the light. The sphere reaches 'radius' units closer and further.
    float zNear = zPullback - radius;
    float zFar  = zPullback + radius;

    // To catch tall buildings or trees BEHIND the camera (outside the frustum), we drastically pull the near plane back.
    // Making zNear negative is perfectly fine in OpenGL ortho projections.
    zNear -= 200.0f; 

    Matrix4x4 lightProjection = Matrix4x4::ortho(minX, maxX, minY, maxY, zNear, zFar);

    return lightProjection * lightView;
}


} // namespace Lizeral
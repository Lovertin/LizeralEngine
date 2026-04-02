// #include "RenderingSystem.h"
// #include "runtime/function/ecs/components/componentAll.h"
// #include "runtime/function/res_type/Material/PBRMaterial.h"
// #include "runtime/resource/resourceManager/resourceManager.h"
// #include "runtime/core/math/matrix4.h"
// #include "runtime/core/math/vector3.h"
// // #include "runtime/function/physics/physicsDebug/PhysicsDebugDrawer.h"

// namespace Lizeral {

//     void RenderingSystem::Initialize() {

//         m_shadowShader = std::make_shared<Shader>(
//             "C:\\Lizeral Engine\\LizeralEngine0.0.1\\engine\\source\\runtime\\function\\sandbox\\shaderTest\\shadow.vert", 
//             "C:\\Lizeral Engine\\LizeralEngine0.0.1\\engine\\source\\runtime\\function\\sandbox\\shaderTest\\shadow.frag"
//         );

//         m_skyboxShader = std::make_shared<Shader>(
//             "C:\\Lizeral Engine\\LizeralEngine0.0.1\\engine\\source\\runtime\\function\\sandbox\\shaderTest\\skybox.vert", 
//             "C:\\Lizeral Engine\\LizeralEngine0.0.1\\engine\\source\\runtime\\function\\sandbox\\shaderTest\\skybox.frag"
//         );

//         SetupSkyboxCube();

//         m_skyboxTexture = ResourceManager::GetInstance().Load<Texture2D>("C:\\Lizeral Engine\\LizeralEngine0.0.1\\asset\\citrus_orchard_road_puresky_4k.hdr");

//         glGenFramebuffers(1, &m_depthMapFBO);

//         glGenTextures(1, &m_depthMapArray);

//         glBindTexture(GL_TEXTURE_2D_ARRAY, m_depthMapArray); 

//         glTexImage3D(GL_TEXTURE_2D_ARRAY, 0, GL_DEPTH_COMPONENT32F, 
//                      SHADOW_WIDTH, SHADOW_HEIGHT, NUM_CASCADES, 
//                      0, GL_DEPTH_COMPONENT, GL_FLOAT, NULL);

//         glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
//         glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
//         glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
//         glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
//         float borderColor[] = { 1.0f, 1.0f, 1.0f, 1.0f };
//         glTexParameterfv(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_BORDER_COLOR, borderColor);

//         glBindFramebuffer(GL_FRAMEBUFFER, m_depthMapFBO);

//         glFramebufferTexture(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, m_depthMapArray, 0);

//         glDrawBuffer(GL_NONE);
//         glReadBuffer(GL_NONE);

//         if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
//             std::cerr << "[RenderingSystem] Error: Shadow Map Framebuffer is not complete!" << std::endl;
//         } else {
//             std::cout << "[RenderingSystem] CSM FBO Initialized Successfully. (Array Layers: " << NUM_CASCADES << ")" << std::endl;
//         }

//         glBindFramebuffer(GL_FRAMEBUFFER, m_DefaultFBO);
//     }

//     void RenderingSystem::Tick(Registry& registry, float deltaTime) {
//         Matrix4x4 viewMatrix = Matrix4x4::IDENTITY;
//         Matrix4x4 projMatrix = Matrix4x4::IDENTITY;
//         Vector3 camPos(0, 0, 0);
//         const CameraComponent* mainCamPtr = nullptr;

//         auto cameraView = registry.view<CameraComponent, TransformComponent>();
//         for (auto entity : cameraView) {
//             auto& cam = cameraView.get<CameraComponent>(entity);
//             if (cam.isMain()) {
//                 mainCamPtr = &cam;
//                 auto& t = cameraView.get<TransformComponent>(entity);
//                 viewMatrix = cam.getViewMatrix();
//                 projMatrix = cam.getProjectionMatrix();
//                 camPos = t.getPosition();
//                 break; 
//             }
//         }

//         if (!mainCamPtr) return; 

//         Vector3 lightDirToSun(1.0f, 1.5f, 1.0f);
//         lightDirToSun.normalise();
//         Vector3 lightColor(3.0f, 3.0f, 3.0f); 
//         auto sunlightView = registry.view<TransformComponent, DirectionLightComponent>();
//         for (auto entity : sunlightView) {
//             auto& light = sunlightView.get<DirectionLightComponent>(entity);
//             if (light.isGlobal()) {
//                 auto& sunTransform = registry.get<TransformComponent>(entity);
//                 lightDirToSun = sunTransform.getForward(); 

//                 lightColor = light.getColor() * light.getIntensity(); 
//                 break;
//             }
//         }

//         std::vector<float> shadowCascadeLevels{ 15.0f, 50.0f, 150.0f }; 
//         std::vector<Matrix4x4> lightSpaceMatrices;
//         lightSpaceMatrices.reserve(4);

//         lightSpaceMatrices.push_back(getLightSpaceMatrix(mainCamPtr->getzNear(), shadowCascadeLevels[0], *mainCamPtr, viewMatrix, lightDirToSun));
//         lightSpaceMatrices.push_back(getLightSpaceMatrix(shadowCascadeLevels[0], shadowCascadeLevels[1], *mainCamPtr, viewMatrix, lightDirToSun));
//         lightSpaceMatrices.push_back(getLightSpaceMatrix(shadowCascadeLevels[1], shadowCascadeLevels[2], *mainCamPtr, viewMatrix, lightDirToSun));
//         lightSpaceMatrices.push_back(getLightSpaceMatrix(shadowCascadeLevels[2], mainCamPtr->getzFar(), *mainCamPtr, viewMatrix, lightDirToSun));

//         glBindFramebuffer(GL_FRAMEBUFFER, m_depthMapFBO);
//         glViewport(0, 0, SHADOW_WIDTH, SHADOW_HEIGHT); 

//         glEnable(GL_DEPTH_TEST); 
//         glDepthMask(GL_TRUE);
//         glEnable(GL_CULL_FACE);
//         glCullFace(GL_BACK); 

//         m_shadowShader->Bind();

//         auto modelView = registry.view<TransformComponent, ModelComponent>();

//         for (int layer = 0; layer < lightSpaceMatrices.size(); ++layer) {
//             glFramebufferTextureLayer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, m_depthMapArray, 0, layer);
//             glClear(GL_DEPTH_BUFFER_BIT); 

//             m_shadowShader->SetUniformMat4f("lightSpaceMatrix", lightSpaceMatrices[layer]);

//             for (auto entity : modelView) {
//                 auto& t = modelView.get<TransformComponent>(entity);
//                 auto& m = modelView.get<ModelComponent>(entity);
//                 if (!m.m_Model) continue;

//                 m_shadowShader->SetUniformMat4f("model", t.getMatrix());
//                 for (auto& mesh : m.m_Model->GetMeshes()) {
//                     mesh.Draw(); 
//                 }
//             }
//         }

//         glBindFramebuffer(GL_FRAMEBUFFER, m_DefaultFBO); 
//         glViewport(m_viewX, m_viewY, m_viewW, m_viewH);

//         glEnable(GL_DEPTH_TEST);
//         glCullFace(GL_BACK); 
        
//         for (auto entity : modelView) {
//             auto& t = modelView.get<TransformComponent>(entity);
//             auto& m = modelView.get<ModelComponent>(entity);

//             if (m.m_ModelPath != m.m_LastLoadedModelPath) {
//                 m.LoadResources(); 
//             }

//             if (!m.m_Model) continue;

//             Matrix4x4 modelMatrix = t.getMatrix();
//             auto& meshes = m.m_Model->GetMeshes();

//             for (size_t i = 0; i < meshes.size(); ++i) {
//                 auto& mesh = meshes[i];

//                 std::shared_ptr<Lizeral::Material> activeMat = mesh.m_Material; 
                
//                 if (i < m.m_OverrideMaterials.size() && m.m_OverrideMaterials[i] != nullptr) {
//                     activeMat = m.m_OverrideMaterials[i]; 
//                 }

//                 auto pbrMat = std::dynamic_pointer_cast<PBRMaterial>(activeMat);
//                 if (!pbrMat || !pbrMat->GetShader()) continue; 

//                 auto shader = pbrMat->GetShader();
//                 shader->Bind();

//                 shader->SetUniformMat4f("model", modelMatrix);
//                 shader->SetUniformMat4f("view", viewMatrix);
//                 shader->SetUniformMat4f("projection", projMatrix);
//                 shader->SetUniformVector3f("u_camPos", camPos);
                
//                 shader->SetUniformVector3f("u_lightDir", lightDirToSun);
//                 shader->SetUniformVector3f("u_lightColor", lightColor);

//                 shader->SetUniform1i("u_cascadeCount", lightSpaceMatrices.size());
//                 for (size_t j = 0; j < lightSpaceMatrices.size(); ++j) {
//                     shader->SetUniformMat4f("u_lightSpaceMatrices[" + std::to_string(j) + "]", lightSpaceMatrices[j]);
//                     if (j < shadowCascadeLevels.size()) {
//                         shader->SetUniform1f("u_cascadePlaneDistances[" + std::to_string(j) + "]", shadowCascadeLevels[j]);
//                     }
//                 }

//                 glActiveTexture(GL_TEXTURE7); 
//                 glBindTexture(GL_TEXTURE_2D_ARRAY, m_depthMapArray); 
//                 shader->SetUniform1i("shadowMap", 7);

//                 pbrMat->BindAndApply();
//                 mesh.Draw(); 
//             }
//         }

//         if (m_cubeVAO != 0 && m_skyboxTexture) {
//             glDepthFunc(GL_LEQUAL); 
//             glDisable(GL_CULL_FACE); 

//             m_skyboxShader->Bind();
//             m_skyboxShader->SetUniformMat4f("view", viewMatrix);
//             m_skyboxShader->SetUniformMat4f("projection", projMatrix);

//             m_skyboxTexture->Bind(0); 
//             m_skyboxShader->SetUniform1i("equirectangularMap", 0);

//             glBindVertexArray(m_cubeVAO);
//             glDrawArrays(GL_TRIANGLES, 0, 36);
//             glBindVertexArray(0);

//             glDepthFunc(GL_LESS); 
//             glEnable(GL_CULL_FACE);
//         }

//         glUseProgram(0);
//     }

//     void RenderingSystem::Shutdown() {
//         if (m_depthMapFBO != 0) {
//             glDeleteFramebuffers(1, &m_depthMapFBO);
//             m_depthMapFBO = 0;
//         }
//     }


//     std::vector<Vector3> getFrustumCornersWorldSpace(const Matrix4x4& projMatrix, const Matrix4x4& viewMatrix) {
//         Matrix4x4 viewProj = projMatrix * viewMatrix;
//         Matrix4x4 invViewProj = viewProj.inverse();

//         std::vector<Vector3> corners;
//         corners.reserve(8);

//         for (unsigned int x = 0; x < 2; ++x) {
//             for (unsigned int y = 0; y < 2; ++y) {
//                 for (unsigned int z = 0; z < 2; ++z) {

//                     float ptX = 2.0f * x - 1.0f;
//                     float ptY = 2.0f * y - 1.0f;
//                     float ptZ = 2.0f * z - 1.0f;

//                     Vector3 ptNDC(ptX, ptY, ptZ);

//                     Vector3 ptWorld = invViewProj * ptNDC; 

//                     Vector4 ptWorld4 = invViewProj * Vector4(ptX, ptY, ptZ, 1.0f);
                    
//                     corners.push_back(Vector3(ptWorld4.x / ptWorld4.w
//                                             , ptWorld4.y / ptWorld4.w, 
//                                               ptWorld4.z / ptWorld4.w));

//                 }
//             }
//         }
//         return corners;
//     }

//     Matrix4x4 RenderingSystem::getLightSpaceMatrix(
//     const float nearPlane, 
//     const float farPlane, 
//     const CameraComponent& cam, 
//     const Matrix4x4& camView, 
//     const Vector3& lightDir) 
// {
//     // 1. USE THE EXACT CAMERA PROJECTION to prevent "backward/mirrored" frustums
//     Matrix4x4 proj = const_cast<CameraComponent&>(cam).buildPerspective(cam.getFov(), cam.getAspect(), nearPlane, farPlane);
//     std::vector<Vector3> corners = getFrustumCornersWorldSpace(proj, camView);

//     // 2. Find the center of the sub-frustum
//     Vector3 center(0.0f, 0.0f, 0.0f);
//     for (const auto& v : corners) {
//         center += v;
//     }
//     center /= corners.size();

//     // 3. Calculate the Bounding Sphere Radius 
//     // This perfectly encapsulates the frustum and stops shadows from "wobbling" when the camera rotates!
//     float radius = 0.0f;
//     for (const auto& v : corners) {
//         float distance = (v - center).length();
//         radius = std::max(radius, distance);
//     }
    
//     // Optional: Snap radius to texel increments here if you want absolutely 0 sub-pixel jitter, 
//     // but the bounding sphere alone fixes 95% of the movement artifacts.

//     // 4. Position the Sun Camera
//     // We pull the sun back by a fixed 'zPullback' distance so it can see objects casting shadows from behind the camera.
//     float zPullback = 300.0f;
//     Vector3 lightPos = center + lightDir * zPullback;
//     Vector3 upVec = (std::abs(lightDir.y) > 0.999f) ? Vector3(0.0f, 0.0f, 1.0f) : Vector3(0.0f, 1.0f, 0.0f);
//     Matrix4x4 lightView = Matrix4x4::lookAt(lightPos, center, upVec);

//     // 5. Build the Ortho Matrix tightly around the sphere
//     // Because the light looks EXACTLY at 'center', the sphere perfectly maps to X/Y bounds from -radius to +radius.
//     // ZERO padding is needed! Your resolution will be incredibly sharp.
//     float minX = -radius;
//     float maxX =  radius;
//     float minY = -radius;
//     float maxY =  radius;

//     // 6. Perfect Z-Bounds
//     // The center is 'zPullback' units away from the light. The sphere reaches 'radius' units closer and further.
//     float zNear = zPullback - radius;
//     float zFar  = zPullback + radius;

//     // To catch tall buildings or trees BEHIND the camera (outside the frustum), we drastically pull the near plane back.
//     // Making zNear negative is perfectly fine in OpenGL ortho projections.
//     zNear -= 200.0f; 

//     Matrix4x4 lightProjection = Matrix4x4::ortho(minX, maxX, minY, maxY, zNear, zFar);

//     return lightProjection * lightView;
// }


// void RenderingSystem::SetupSkyboxCube() {
//     float skyboxVertices[] = {
//         // positions          
//         -1.0f,  1.0f, -1.0f, -1.0f, -1.0f, -1.0f,  1.0f, -1.0f, -1.0f,
//             1.0f, -1.0f, -1.0f,  1.0f,  1.0f, -1.0f, -1.0f,  1.0f, -1.0f,
//         -1.0f, -1.0f,  1.0f, -1.0f, -1.0f, -1.0f, -1.0f,  1.0f, -1.0f,
//         -1.0f,  1.0f, -1.0f, -1.0f,  1.0f,  1.0f, -1.0f, -1.0f,  1.0f,
//             1.0f, -1.0f, -1.0f,  1.0f, -1.0f,  1.0f,  1.0f,  1.0f,  1.0f,
//             1.0f,  1.0f,  1.0f,  1.0f,  1.0f, -1.0f,  1.0f, -1.0f, -1.0f,
//         -1.0f, -1.0f,  1.0f, -1.0f,  1.0f,  1.0f,  1.0f,  1.0f,  1.0f,
//             1.0f,  1.0f,  1.0f,  1.0f, -1.0f,  1.0f, -1.0f, -1.0f,  1.0f,
//         -1.0f,  1.0f, -1.0f,  1.0f,  1.0f, -1.0f,  1.0f,  1.0f,  1.0f,
//             1.0f,  1.0f,  1.0f, -1.0f,  1.0f,  1.0f, -1.0f,  1.0f, -1.0f,
//         -1.0f, -1.0f, -1.0f, -1.0f, -1.0f,  1.0f,  1.0f, -1.0f, -1.0f,
//             1.0f, -1.0f, -1.0f, -1.0f, -1.0f,  1.0f,  1.0f, -1.0f,  1.0f
//     };
//     glGenVertexArrays(1, &m_cubeVAO);
//     glGenBuffers(1, &m_cubeVBO);
//     glBindVertexArray(m_cubeVAO);
//     glBindBuffer(GL_ARRAY_BUFFER, m_cubeVBO);
//     glBufferData(GL_ARRAY_BUFFER, sizeof(skyboxVertices), &skyboxVertices, GL_STATIC_DRAW);
//     glEnableVertexAttribArray(0);
//     glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
//     glBindVertexArray(0);
// }

// } // namespace Lizeral
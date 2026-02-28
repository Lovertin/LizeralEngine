#include <iostream>
#include <vector>
#include <cmath>
#include <memory>

// 1. OpenGL Headers
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <QApplication>

// 2. Engine Headers
#include "runtime/function/ecs/registry.h"
#include "runtime/function/ecs/components/Collider/ColliderComponent.h"
#include "runtime/function/ecs/components/RigidBody/RigidBodyComponent.h"
#include "runtime/function/ecs/components/Transform/TransformComponent.h"
#include "runtime/function/ecs/components/Model/ModelComponent.h"
#include "runtime/function/ecs/components/Light/DirectionalLightComponent.h"
#include "runtime/function/physics/PhysicsScene.h"
#include "runtime/function/physics/PhysicsSystem.h"
#include "runtime/function/render/RenderingSystem/RenderingSystem.h"
#include "runtime/core/math/vector3.h"
#include "runtime/core/math/quaternion.h"
#include "runtime/core/math/matrix4.h"

// 3. 射线与命中数据结构
#include "runtime/function/physics/phsicsEntityheaders.h"

// 4. 相机与资源系统
#include "runtime/function/input/input.h"
#include "runtime/function/ecs/components/Camera/CameraComponent.h"
#include "runtime/function/ecs/components/Camera/CameraControlComponent.h"
#include "runtime/function/render/CameraControlSystem/CameraControlSystem.h"
#include "runtime/function/render/CameraSystem/CameraSystem.h"
#include "runtime/resource/resourceManager/resourceManager.h"
#include "runtime/function/res_type/shader/shader.h"
#include "runtime/function/res_type/Material/Material.h"
#include "runtime/function/res_type/Material/PBRMaterial.h"
#include "runtime/function/res_type/texture/TextureCube.h"
#include "runtime/function/res_type/Model/Mesh.h"
#include "runtime/function/res_type/Model/Model.h"

using namespace Lizeral;

const int WIN_WIDTH = 1280;
const int WIN_HEIGHT = 720;

// --- 辅助：将 Lizeral Matrix4x4 转换为 OpenGL 格式并加载 (用于物理 Debug 线框) ---
void LoadEngineMatrixToOpenGL(GLenum mode, const Matrix4x4& m) {
    float glMat[16];
    glMat[0] = m[0][0]; glMat[1] = m[1][0]; glMat[2] = m[2][0]; glMat[3] = m[3][0];
    glMat[4] = m[0][1]; glMat[5] = m[1][1]; glMat[6] = m[2][1]; glMat[7] = m[3][1];
    glMat[8] = m[0][2]; glMat[9] = m[1][2]; glMat[10]= m[2][2]; glMat[11]= m[3][2];
    glMat[12]= m[0][3]; glMat[13]= m[1][3]; glMat[14]= m[2][3]; glMat[15]= m[3][3];

    glMatrixMode(mode);
    glLoadMatrixf(glMat);
}

void framebuffer_size_callback(GLFWwindow* window, int width, int height) {
    glViewport(0, 0, width, height);
}

int main() {
    // ========================================================
    // 1. 初始化 GLFW、GLAD、系统
    // ========================================================
    if (!glfwInit()) return -1;
    GLFWwindow* window = glfwCreateWindow(WIN_WIDTH, WIN_HEIGHT, "Lizeral Engine - PBR & Shadow", NULL, NULL);
    if (!window) { glfwTerminate(); return -1; }
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) return -1;

    Registry registry;
    
    PhysicsScene physicsScene;
    PhysicsSystem physicsSystem;
    physicsScene.Initialize();
    physicsSystem.Initialize(&physicsScene);
    
    RenderingSystem renderSystem;
    renderSystem.Initialize();

    Input::GetInstance().Init(window);
    CameraSystem cameraSystem;
    CameraControlSystem cameraControlSystem;

    std::cout << "All Systems Initialized. Hold RMB + WASD to move." << std::endl;

    // ========================================================
    // 2. 资源加载与材质配置
    // ========================================================
    ResourceManager::GetInstance().SetRootPath("");
    
    auto shuttleModel = ResourceManager::GetInstance().Load<Model>("C:\\Lizeral Engine\\LizeralEngine0.0.1\\asset\\mazda_glb.glb");
    auto boxModel = ResourceManager::GetInstance().Load<Model>("C:\\Lizeral Engine\\LizeralEngine0.0.1\\asset\\cardboard_box.glb");
    auto pureBox = ResourceManager::GetInstance().Load<Model>("C:\\Lizeral Engine\\LizeralEngine0.0.1\\asset\\box_with_uv.glb");

    auto pbrShader = std::make_shared<Shader>(
        "C:\\Lizeral Engine\\LizeralEngine0.0.1\\engine\\source\\runtime\\function\\sandbox\\shaderTest\\solid.vert", 
        "C:\\Lizeral Engine\\LizeralEngine0.0.1\\engine\\source\\runtime\\function\\sandbox\\shaderTest\\solid.frag"
    );
    auto skyboxShader = std::make_shared<Shader>(
        "C:\\Lizeral Engine\\LizeralEngine0.0.1\\engine\\source\\runtime\\function\\sandbox\\shaderTest\\skybox.vert",
        "C:\\Lizeral Engine\\LizeralEngine0.0.1\\engine\\source\\runtime\\function\\sandbox\\shaderTest\\skybox.frag"
    );

    auto irradianceMap = ResourceManager::GetInstance().Load<TextureCube>("C:\\Lizeral Engine\\LizeralEngine0.0.1\\asset\\texture\\skyCube\\skybox_irradiance");
    auto specularMap = ResourceManager::GetInstance().Load<TextureCube>("C:\\Lizeral Engine\\LizeralEngine0.0.1\\asset\\texture\\skyCube\\skybox_specular");

    // 统一配置马自达的 PBR 材质
    for (auto& mesh : shuttleModel->GetMeshes()) {
        auto pbrMat = std::dynamic_pointer_cast<PBRMaterial>(mesh.m_Material);
        if (pbrMat) {
            pbrMat->SetShader(pbrShader);
            pbrMat->m_IrradianceMap = irradianceMap;
            pbrMat->m_PrefilterMap = specularMap;
            pbrMat->m_Metallic = 0.2f;   
            pbrMat->m_Roughness = 0.15f; 
        }
    }

    // 统一配置盒子的 PBR 材质 (包括地面、下落方块、头顶悬浮方块都会共用这个属性)
    for (auto& mesh : boxModel->GetMeshes()) {
        auto pbrMat = std::dynamic_pointer_cast<PBRMaterial>(mesh.m_Material);
        if (pbrMat) {
            pbrMat->SetShader(pbrShader);
            // 这里故意注释掉 IBL，让盒子和地面呈现纯粹的泥土亚光感，不受天空盒反射干扰
            // pbrMat->m_IrradianceMap = irradianceMap;
            // pbrMat->m_PrefilterMap = specularMap;
            pbrMat->m_Metallic = 0.0f;   
            pbrMat->m_Roughness = 0.9f; 
        }
    }

    for (auto& mesh : pureBox->GetMeshes()) {
        auto pbrMat = std::dynamic_pointer_cast<PBRMaterial>(mesh.m_Material);
        if (pbrMat) {
            pbrMat->SetShader(pbrShader);
            // 这里故意注释掉 IBL，让盒子和地面呈现纯粹的泥土亚光感，不受天空盒反射干扰
            // pbrMat->m_IrradianceMap = irradianceMap;
            // pbrMat->m_PrefilterMap = specularMap;
            pbrMat->m_Metallic = 0.1f;   
            pbrMat->m_Roughness = 0.9f; 
        }
    }

    // ========================================================
    // 3. 构建 ECS 场景实体
    // ========================================================

    // --- A. 主相机 ---
    Entity cameraEntity = registry.create();
    {
        auto& t = registry.emplace<TransformComponent>(cameraEntity);
        t.setPosition(Vector3(0.0f, 10.0f, 20.0f)); 
        
        auto& cam = registry.emplace<CameraComponent>(cameraEntity);
        cam.setFov(45.0f);
        cam.setAspect((float)WIN_WIDTH / (float)WIN_HEIGHT);
        cam.setzNear(0.1f);
        cam.setzFar(1000.0f);
        cam.setMain(true);

        auto& ctrl = registry.emplace<CameraControlComponent>(cameraEntity);
        ctrl.setMoveSpeed(2.0f);
        ctrl.setSensitivityX(0.1f); 
        ctrl.setSensitivityY(0.1f);
        ctrl.setYaw(-90.0f);
    }

    // // --- B. 地面 (升级为 PBR Box 模型！) ---
    Entity ground = registry.create();
    {
        auto& t = registry.emplace<TransformComponent>(ground);
        t.setPosition(Vector3(0, 0, 0));
        t.setScale(Vector3(30.0f, 1.0f, 30.0f)); // 压扁并放大
        
        // 【关键】：地面现在也是真正的模型了，能够接收阴影！
        registry.emplace<ModelComponent>(ground, pureBox);

        auto& c = registry.emplace<ColliderComponent>(ground);
        c.setType(ColliderType::Box);
        c.setSize(Vector3(2.0f, 2.0f, 2.0f));
        auto& rb = registry.emplace<RigidBodyComponent>(ground);
        rb.setMass(0.0f); rb.setFriction(0.8f); rb.setRestitution(1.0f);
    }

    // --- C. 动态方块 (同样升级为 PBR 模型) ---
    for (int y = 0; y < 4; y++) {
        for (int x = 0; x < 3; x++) {
            Entity cube = registry.create();
            auto& t = registry.emplace<TransformComponent>(cube);
            t.setPosition(Vector3(x * 2.5f - 2.5f, 10.0f + y * 3.0f, 0)); 
            t.setScale(Vector3(1, 1, 1));
            
            registry.emplace<ModelComponent>(cube, pureBox);

            auto& c = registry.emplace<ColliderComponent>(cube);
            c.setType(ColliderType::Box);
            c.setSize(Vector3(2.0f, 2.0f, 2.0f));
            auto& rb = registry.emplace<RigidBodyComponent>(cube);
            rb.setMass(1.0f); rb.setFriction(0.5f); rb.setRestitution(0.7f);
        }
    }

    // --- D. 纯渲染锚点：马自达 (无物理组件) ---
    Entity carEntity = registry.create();
    {
        auto& t = registry.emplace<TransformComponent>(carEntity);
        t.setPosition(Vector3(20.0f, 10.0f, 0.0f)); 
        t.setScale(Vector3(2.0f, 2.0f, 2.0f));
        registry.emplace<ModelComponent>(carEntity, shuttleModel);

        // auto& rb = registry.emplace<RigidBodyComponent>(carEntity);
        // rb.setMass(1.0f); rb.setFriction(0.5f); rb.setRestitution(0.5f);

        // auto& c = registry.emplace<ColliderComponent>(carEntity);
        // c.setType(ColliderType::Box);
        // c.setSize(Vector3(2.0f, 2.0f, 2.0f));
    }

    // --- E. 纯渲染锚点：悬浮巨盒 (无物理组件) ---
    Entity floatingBoxEntity = registry.create();
    {
        auto& t = registry.emplace<TransformComponent>(floatingBoxEntity);
        // 放在车的正上方
        t.setPosition(Vector3(20.0f, 7.0f, 0.0f)); 
        t.setScale(Vector3(0.1f, 0.1f, 0.1f)); 
        registry.emplace<ModelComponent>(floatingBoxEntity, boxModel);

        auto& rb = registry.emplace<RigidBodyComponent>(floatingBoxEntity);
        rb.setMass(1.0f); rb.setFriction(0.5f); rb.setRestitution(0.5f);

        auto& c = registry.emplace<ColliderComponent>(floatingBoxEntity);
        c.setType(ColliderType::Box);
        c.setSize(Vector3(80.0f, 40.0f, 80.0f));
    }
    // --- F. 定向光太阳 ---
    Entity SunLight =registry.create();
    {
        auto& t=registry.emplace<TransformComponent>(SunLight);
        t.setForward(Vector3(1.0f, 1.0f, 1.0f));

        auto& l =registry.emplace<DirectionLightComponent>(SunLight);
        l.setIntensity(3.0);
    }

    // --- 天空盒 VAO 初始化 ---
    float skyboxVertices[] = {
        -1.0f,  1.0f, -1.0f, -1.0f, -1.0f, -1.0f,  1.0f, -1.0f, -1.0f,  1.0f, -1.0f, -1.0f,  1.0f,  1.0f, -1.0f, -1.0f,  1.0f, -1.0f,
        -1.0f, -1.0f,  1.0f, -1.0f, -1.0f, -1.0f, -1.0f,  1.0f, -1.0f, -1.0f,  1.0f, -1.0f, -1.0f,  1.0f,  1.0f, -1.0f, -1.0f,  1.0f,
         1.0f, -1.0f, -1.0f,  1.0f, -1.0f,  1.0f,  1.0f,  1.0f,  1.0f,  1.0f,  1.0f,  1.0f,  1.0f,  1.0f, -1.0f,  1.0f, -1.0f, -1.0f,
        -1.0f, -1.0f,  1.0f, -1.0f,  1.0f,  1.0f,  1.0f,  1.0f,  1.0f,  1.0f,  1.0f,  1.0f,  1.0f, -1.0f,  1.0f, -1.0f, -1.0f,  1.0f,
        -1.0f,  1.0f, -1.0f,  1.0f,  1.0f, -1.0f,  1.0f,  1.0f,  1.0f,  1.0f,  1.0f,  1.0f, -1.0f,  1.0f,  1.0f, -1.0f,  1.0f, -1.0f,
        -1.0f, -1.0f, -1.0f, -1.0f, -1.0f,  1.0f,  1.0f, -1.0f, -1.0f,  1.0f, -1.0f, -1.0f, -1.0f, -1.0f,  1.0f,  1.0f, -1.0f,  1.0f
    };
    unsigned int skyboxVAO, skyboxVBO;
    glGenVertexArrays(1, &skyboxVAO);
    glGenBuffers(1, &skyboxVBO);
    glBindVertexArray(skyboxVAO);
    glBindBuffer(GL_ARRAY_BUFFER, skyboxVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(skyboxVertices), &skyboxVertices, GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glBindVertexArray(0);

    // ========================================================
    // 4. 主循环
    // ========================================================
    float lastTime = 0.0f;
    glEnable(GL_DEPTH_TEST);

    while (!glfwWindowShouldClose(window)) {
        float currentTime = (float)glfwGetTime();
        float deltaTime = currentTime - lastTime;
        lastTime = currentTime;

        // --- A. 输入与 Raycast 逻辑 ---
        Input::GetInstance().Tick();
        if (Input::GetInstance().GetKey(Key::ESC)) glfwSetWindowShouldClose(window, true);

        if (Input::GetInstance().GetMouseButtonDown(MouseButton::Left)) {
            auto& mainCam = registry.get<CameraComponent>(cameraEntity);
            auto& camTrans = registry.get<TransformComponent>(cameraEntity);
            Vector2 mousePos = Input::GetInstance().GetMousePosition();
            float ndcX = (2.0f * mousePos.x) / WIN_WIDTH - 1.0f;
            float ndcY = 1.0f - (2.0f * mousePos.y) / WIN_HEIGHT;
            Matrix4x4 invProj = mainCam.getProjectionMatrix().inverse();
            Matrix4x4 invView = mainCam.getViewMatrix().inverse();
            Vector3 nearPointWorld = invView * (invProj * Vector3(ndcX, ndcY, -1.0f));
            Vector3 farPointWorld = invView * (invProj * Vector3(ndcX, ndcY, 1.0f));
            Vector3 rayDir = (farPointWorld - nearPointWorld).normalisedCopy();
            
            RaycastHit hitInfo;
            if (physicsSystem.Raycast(Ray(camTrans.getPosition(), rayDir), hitInfo)) {
                std::cout << "[Hit] ID: " << hitInfo.entity << " Pos: " << hitInfo.point.x << ", " << hitInfo.point.y << ", " << hitInfo.point.z << "\n";
            }
        }

        // --- B. 物理与系统模拟 ---
        physicsSystem.Tick(deltaTime, registry);
        cameraControlSystem.Tick(registry, window); 
        cameraSystem.Tick(registry);

        // ====================================================
        // --- C. 核心渲染管线启动 ---
        // ====================================================
        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // 1. PBR 与 阴影渲染系统接管 (全自动遍历所有 ModelComponent)
        renderSystem.Tick(registry, deltaTime);

        // 2. 补画天空盒
        // 【防污染】：切回 0 号纹理，防止 PBR 管线的深度图残留在纹理槽里
        glActiveTexture(GL_TEXTURE0); 
        glDepthFunc(GL_LEQUAL); 
        skyboxShader->Bind();

        auto& mainCam = registry.get<CameraComponent>(cameraEntity);
        Matrix4x4 viewMatrix = mainCam.getViewMatrix();
        Matrix4x4 viewNoTranslation = Matrix4x4::IDENTITY;
        for(int i=0; i<3; ++i) for(int j=0; j<3; ++j) viewNoTranslation[i][j] = viewMatrix[i][j];

        skyboxShader->SetUniformMat4f("view", viewNoTranslation);
        skyboxShader->SetUniformMat4f("projection", mainCam.getProjectionMatrix());

        glBindVertexArray(skyboxVAO);
        specularMap->Bind(0); 
        skyboxShader->SetUniform1i("skybox", 0);
        glDrawArrays(GL_TRIANGLES, 0, 36);
        
        glBindVertexArray(0);
        glDepthFunc(GL_LESS); 
        glUseProgram(0);

        // 3. 物理 Debug X光机 (绘制红色/绿色/白色的线框边界)
        glDisable(GL_CULL_FACE); // 确保各种刁钻角度都能看到物理线框
        LoadEngineMatrixToOpenGL(GL_PROJECTION, mainCam.getProjectionMatrix());
        LoadEngineMatrixToOpenGL(GL_MODELVIEW, mainCam.getViewMatrix());
        
        // physicsSystem.DebugDrawWorld();
        
        glMatrixMode(GL_PROJECTION); glLoadIdentity();
        glMatrixMode(GL_MODELVIEW); glLoadIdentity();

        glfwSwapBuffers(window);
    }

    // 5. 清理
    physicsSystem.Shutdown();
    renderSystem.Shutdown();
    glfwTerminate();

    return 0;
}
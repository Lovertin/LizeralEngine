#include <iostream>
#include <vector>
#include <cmath>

// 1. OpenGL Headers
#include <glad/glad.h>
#include <GLFW/glfw3.h>

// 2. Engine Headers
#include "runtime/function/ecs/registry.h"
#include "runtime/function/ecs/components/Collider/ColliderComponent.h"
#include "runtime/function/ecs/components/RigidBody/RigidBodyComponent.h"
#include "runtime/function/ecs/components/Transform/TransformComponent.h"
#include "runtime/function/physics/PhysicsScene.h"
#include "runtime/function/physics/PhysicsSystem.h"
#include "runtime/core/math/vector3.h"
#include "runtime/core/math/quaternion.h"
#include "runtime/core/math/matrix4.h"

// --- 新增：相机与输入系统头文件 ---
#include "runtime/function/input/input.h"
#include "runtime/function/ecs/components/Camera/CameraComponent.h"
#include "runtime/function/ecs/components/Camera/CameraControlComponent.h"
#include "runtime/function/render/CameraControlSystem/CameraControlSystem.h"
#include "runtime/function/render/CameraSystem/CameraSystem.h"

using namespace Lizeral;

// --- 全局变量 ---
const int WIN_WIDTH = 1280;
const int WIN_HEIGHT = 720;

// --- 辅助：将 Lizeral Matrix4x4 转换为 OpenGL 格式并加载 ---
// OpenGL 使用列主序 (Column-Major)，C++ 数学库通常是行主序 (Row-Major)
// 如果你的矩阵类内存布局是行主序，这里需要转置
void LoadEngineMatrixToOpenGL(GLenum mode, const Matrix4x4& m) {
    float glMat[16];

    // 假设 Matrix4x4 提供了按 (row, col) 访问或 flat 数组
    // 这里进行手动转置以确保万无一失： m[row][col] -> glMat[col * 4 + row]
    
    // Column 0
    glMat[0] = m[0][0]; glMat[1] = m[1][0]; glMat[2] = m[2][0]; glMat[3] = m[3][0];
    // Column 1
    glMat[4] = m[0][1]; glMat[5] = m[1][1]; glMat[6] = m[2][1]; glMat[7] = m[3][1];
    // Column 2
    glMat[8] = m[0][2]; glMat[9] = m[1][2]; glMat[10]= m[2][2]; glMat[11]= m[3][2];
    // Column 3
    glMat[12]= m[0][3]; glMat[13]= m[1][3]; glMat[14]= m[2][3]; glMat[15]= m[3][3];

    glMatrixMode(mode);
    glLoadMatrixf(glMat);
}

// --- 辅助：将 Lizeral Transform 转换为 OpenGL 矩阵 ---
void ApplyTransform(const TransformComponent& t) {
    glTranslatef(t.getPosition().x, t.getPosition().y, t.getPosition().z);
    
    Quaternion q = t.getRotation();
    float xx = q.x * q.x; float yy = q.y * q.y; float zz = q.z * q.z;
    float xy = q.x * q.y; float xz = q.x * q.z; float yz = q.y * q.z;
    float wx = q.w * q.x; float wy = q.w * q.y; float wz = q.w * q.z;

    float mat[16];
    mat[0] = 1.0f - 2.0f * (yy + zz); mat[1] = 2.0f * (xy + wz); mat[2] = 2.0f * (xz - wy); mat[3] = 0.0f;
    mat[4] = 2.0f * (xy - wz); mat[5] = 1.0f - 2.0f * (xx + zz); mat[6] = 2.0f * (yz + wx); mat[7] = 0.0f;
    mat[8] = 2.0f * (xz + wy); mat[9] = 2.0f * (yz - wx); mat[10] = 1.0f - 2.0f * (xx + yy); mat[11] = 0.0f;
    mat[12] = 0.0f; mat[13] = 0.0f; mat[14] = 0.0f; mat[15] = 1.0f;

    glMultMatrixf(mat);
    glScalef(t.getScale().x, t.getScale().y, t.getScale().z);
}

// --- 辅助：画一个简单的立方体 ---
void DrawCube(const Vector3& color, const Vector3& size) {
    glColor3f(color.x, color.y, color.z);
    float hx = size.x * 0.5f;
    float hy = size.y * 0.5f;
    float hz = size.z * 0.5f;

    glBegin(GL_QUADS);
    // ... (保持原有的顶点绘制代码不变) ...
    // Front face
    glVertex3f(-hx, -hy,  hz); glVertex3f( hx, -hy,  hz); glVertex3f( hx,  hy,  hz); glVertex3f(-hx,  hy,  hz);
    // Back face
    glVertex3f(-hx, -hy, -hz); glVertex3f(-hx,  hy, -hz); glVertex3f( hx,  hy, -hz); glVertex3f( hx, -hy, -hz);
    // Top face
    glVertex3f(-hx,  hy, -hz); glVertex3f(-hx,  hy,  hz); glVertex3f( hx,  hy,  hz); glVertex3f( hx,  hy, -hz);
    // Bottom face
    glVertex3f(-hx, -hy, -hz); glVertex3f( hx, -hy, -hz); glVertex3f( hx, -hy,  hz); glVertex3f(-hx, -hy,  hz);
    // Right face
    glVertex3f( hx, -hy, -hz); glVertex3f( hx,  hy, -hz); glVertex3f( hx,  hy,  hz); glVertex3f( hx, -hy,  hz);
    // Left face
    glVertex3f(-hx, -hy, -hz); glVertex3f(-hx, -hy,  hz); glVertex3f(-hx,  hy,  hz); glVertex3f(-hx,  hy, -hz);
    glEnd();
}

// 窗口 Resize 回调
void framebuffer_size_callback(GLFWwindow* window, int width, int height) {
    glViewport(0, 0, width, height);
}

int main() {
    // --------------------------------------------------------
    // 1. 初始化 GLFW 和 GLAD
    // --------------------------------------------------------
    if (!glfwInit()) return -1;
    
    GLFWwindow* window = glfwCreateWindow(WIN_WIDTH, WIN_HEIGHT, "Lizeral Engine - Roaming Camera", NULL, NULL);
    if (!window) { glfwTerminate(); return -1; }
    
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    // --------------------------------------------------------
    // 2. 初始化 ECS 和 系统
    // --------------------------------------------------------
    Registry registry;
    
    // A. 物理系统
    PhysicsScene physicsScene;
    PhysicsSystem physicsSystem;
    physicsScene.Initialize();
    physicsSystem.Initialize(&physicsScene);

    // B. 输入系统 (这一步会隐藏光标并锁定)
    Input::GetInstance().Init(window);

    // C. 相机系统
    CameraSystem cameraSystem;
    CameraControlSystem cameraControlSystem;

    std::cout << "All Systems Initialized. Press WASD to move, Mouse to look, ESC to exit." << std::endl;

    // --------------------------------------------------------
    // 3. 创建场景实体
    // --------------------------------------------------------

    // --- A. 创建主相机 (New!) ---
    Entity cameraEntity = registry.create();
    {
        // 1. Transform: 放在 (0, 10, 20) 
        auto& t = registry.emplace<TransformComponent>(cameraEntity);
        t.setPosition(Vector3(0.0f, 10.0f, 20.0f)); 
        // 初始朝向：默认看向 -Z，所以这里是看向屏幕内。如果不放心，可以先设个初始 Pitch/Yaw
        
        // 2. Camera: 光学属性
        auto& cam = registry.emplace<CameraComponent>(cameraEntity);
        cam.setFov(45.0f);
        cam.setAspect((float)WIN_WIDTH / (float)WIN_HEIGHT);
        cam.setzNear(0.1f);
        cam.setzFar(1000.0f);
        cam.setMain(true);

        // 3. Control: 漫游控制
        auto& ctrl = registry.emplace<CameraControlComponent>(cameraEntity);
        ctrl.setMoveSpeed(1.0f);     // 移动速度
        ctrl.setSensitivityX(0.02f);   // 鼠标灵敏度
        ctrl.setSensitivityY(0.1f);
        ctrl.setYaw(-90.0f);          // 初始 Yaw (-90 指向 -Z)
    }

    // --- B. 地面 ---
    {
        Entity ground = registry.create();
        auto& t = registry.emplace<TransformComponent>(ground);
        t.setPosition(Vector3(0, -10, 0));
        t.setScale(Vector3(1, 1, 1));
        
        auto& c = registry.emplace<ColliderComponent>(ground);
        c.setType(ColliderType::Box);
        c.setSize(Vector3(100.0f, 2.0f, 100.0f)); // 加大一点地面

        auto& rb = registry.emplace<RigidBodyComponent>(ground);
        rb.setMass(0.0f); 
        rb.setFriction(0.8f);
        rb.setRestitution(0.5f);
    }

    // --- C. 动态方块 ---
    for (int y = 0; y < 4; y++) {
        for (int x = 0; x < 3; x++) {
            Entity cube = registry.create();
            auto& t = registry.emplace<TransformComponent>(cube);
            t.setPosition(Vector3(x * 2.5f - 2.5f, 10.0f + y * 3.0f, 0)); 
            t.setScale(Vector3(1, 1, 1));

            auto& c = registry.emplace<ColliderComponent>(cube);
            c.setType(ColliderType::Box);
            c.setSize(Vector3(1.0f,1.0f,1.0f));

            auto& rb = registry.emplace<RigidBodyComponent>(cube);
            rb.setMass(1.0f);
            rb.setFriction(0.5f);
            rb.setRestitution(0.2f);
        }
    }

    // --------------------------------------------------------
    // 4. 主循环
    // --------------------------------------------------------
    float lastTime = 0.0f;
    glEnable(GL_DEPTH_TEST);

    while (!glfwWindowShouldClose(window)) {
        float currentTime = (float)glfwGetTime();
        float deltaTime = currentTime - lastTime;
        lastTime = currentTime;

        // --- 逻辑更新 ---
        
        // 1. 输入处理
        Input::GetInstance().Tick();
        if (Input::GetInstance().GetKey(Key::ESC)) glfwSetWindowShouldClose(window, true);

        // 2. 物理模拟
        physicsSystem.Tick(deltaTime, registry);

        // 3. 相机控制 (读取 Input -> 修改 Transform)
        // 你的 CameraControlSystem 可能需要 delta time 来计算速度，建议加上
        // cameraControlSystem.Tick(registry, deltaTime); 
        // 暂时使用你目前的无参版本，但在内部记得处理 dt
        cameraControlSystem.Tick(registry); 

        // 4. 相机矩阵计算 (读取 Transform -> 更新 CameraComponent 矩阵)
        cameraSystem.Tick(registry);


        // --- 渲染 ---
        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // --- 应用相机矩阵 (核心修改点) ---
        // 我们不再使用 gluPerspective 或 glTranslate/Rotate 手动推导
        // 而是直接从 CameraComponent 拿算好的矩阵
        
        // 为了方便，这里直接获取刚才创建的 cameraEntity 的组件
        // 在完整的引擎中，RenderSystem 会去查找 isMain() 的相机
        auto& mainCam = registry.get<CameraComponent>(cameraEntity);

        // 1. 设置投影矩阵 (Projection)
        LoadEngineMatrixToOpenGL(GL_PROJECTION, mainCam.getProjectionMatrix());

        // 2. 设置视图矩阵 (View)
        LoadEngineMatrixToOpenGL(GL_MODELVIEW, mainCam.getViewMatrix());

        // --- 绘制物体 ---
        auto view = registry.view<TransformComponent, ColliderComponent>();
        
        for (auto entity : view) {
            auto& t = view.get<TransformComponent>(entity);
            auto& c = view.get<ColliderComponent>(entity);

            glPushMatrix();
            
            // 应用物体变换
            ApplyTransform(t);
            glTranslatef(c.getOffset().x, c.getOffset().y, c.getOffset().z);

            // 绘制
            if (c.getType() == ColliderType::Box) {
                if (c.getSize().x > 10.0f) 
                    DrawCube(Vector3(0.2f, 0.8f, 0.2f), c.getSize());
                else 
                    DrawCube(Vector3(1.0f, 0.5f, 0.2f), c.getSize());
            }
            
            glPopMatrix();
        }

        glfwSwapBuffers(window);
        // glfwPollEvents() 已经被 Input::Tick() 内部调用了，这里可以不需要
    }

    // 5. 清理
    physicsSystem.Shutdown();
    glfwTerminate();

    return 0;
}
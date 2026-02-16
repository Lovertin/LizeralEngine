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
#include "runtime/function/physics/PhysicsSystem.h"
#include "runtime/function/physics/PhysicsScene.h"
#include "runtime/core/math/vector3.h"
#include "runtime/core/math/quaternion.h"

using namespace Lizeral;

// --- 全局变量 ---
const int WIN_WIDTH = 1280;
const int WIN_HEIGHT = 720;

// --- 辅助：将 Lizeral Transform 转换为 OpenGL 矩阵 ---
void ApplyTransform(const TransformComponent& t) {
    // 1. 平移
    glTranslatef(t.getPosition().x, t.getPosition().y, t.getPosition().z);

    // 2. 旋转 (Quaternion -> Axis Angle)
    // 这里为了简单，我们手动解算一下四元数到旋转矩阵，或者直接推导
    // OpenGL 的 glMultMatrixf 需要列主序矩阵
    Quaternion q = t.getRotation();
    
    // 标准的四元数转旋转矩阵公式
    float xx = q.x * q.x; float yy = q.y * q.y; float zz = q.z * q.z;
    float xy = q.x * q.y; float xz = q.x * q.z; float yz = q.y * q.z;
    float wx = q.w * q.x; float wy = q.w * q.y; float wz = q.w * q.z;

    float mat[16];
    mat[0] = 1.0f - 2.0f * (yy + zz);
    mat[1] = 2.0f * (xy + wz);
    mat[2] = 2.0f * (xz - wy);
    mat[3] = 0.0f;

    mat[4] = 2.0f * (xy - wz);
    mat[5] = 1.0f - 2.0f * (xx + zz);
    mat[6] = 2.0f * (yz + wx);
    mat[7] = 0.0f;

    mat[8] = 2.0f * (xz + wy);
    mat[9] = 2.0f * (yz - wx);
    mat[10] = 1.0f - 2.0f * (xx + yy);
    mat[11] = 0.0f;

    mat[12] = 0.0f; mat[13] = 0.0f; mat[14] = 0.0f; mat[15] = 1.0f;

    glMultMatrixf(mat);

    // 3. 缩放
    glScalef(t.getScale().x, t.getScale().y, t.getScale().z);
}

// --- 辅助：画一个简单的立方体 ---
void DrawCube(const Vector3& color, const Vector3& size) {
    glColor3f(color.x, color.y, color.z);
    
    // 为了对应 Bullet 的 BoxShape，我们需要画一个中心在原点的盒子
    // Size 是全长，所以顶点是 +/- size/2
    float hx = size.x * 0.5f;
    float hy = size.y * 0.5f;
    float hz = size.z * 0.5f;

    glBegin(GL_QUADS);
    // Front face
    glVertex3f(-hx, -hy,  hz);
    glVertex3f( hx, -hy,  hz);
    glVertex3f( hx,  hy,  hz);
    glVertex3f(-hx,  hy,  hz);
    // Back face
    glVertex3f(-hx, -hy, -hz);
    glVertex3f(-hx,  hy, -hz);
    glVertex3f( hx,  hy, -hz);
    glVertex3f( hx, -hy, -hz);
    // Top face
    glVertex3f(-hx,  hy, -hz);
    glVertex3f(-hx,  hy,  hz);
    glVertex3f( hx,  hy,  hz);
    glVertex3f( hx,  hy, -hz);
    // Bottom face
    glVertex3f(-hx, -hy, -hz);
    glVertex3f( hx, -hy, -hz);
    glVertex3f( hx, -hy,  hz);
    glVertex3f(-hx, -hy,  hz);
    // Right face
    glVertex3f( hx, -hy, -hz);
    glVertex3f( hx,  hy, -hz);
    glVertex3f( hx,  hy,  hz);
    glVertex3f( hx, -hy,  hz);
    // Left face
    glVertex3f(-hx, -hy, -hz);
    glVertex3f(-hx, -hy,  hz);
    glVertex3f(-hx,  hy,  hz);
    glVertex3f(-hx,  hy, -hz);
    glEnd();
}

int main() {
    // --------------------------------------------------------
    // 1. 初始化 GLFW 和 GLAD
    // --------------------------------------------------------
    if (!glfwInit()) return -1;
    
    GLFWwindow* window = glfwCreateWindow(WIN_WIDTH, WIN_HEIGHT, "Lizeral Physics Sandbox", NULL, NULL);
    if (!window) { glfwTerminate(); return -1; }
    
    glfwMakeContextCurrent(window);
    
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    // OpenGL 设置
    glEnable(GL_DEPTH_TEST);
    glViewport(0, 0, WIN_WIDTH, WIN_HEIGHT);

    // --------------------------------------------------------
    // 2. 初始化 ECS 和 物理系统
    // --------------------------------------------------------
    Registry registry;
    PhysicsScene physicsScene;
    PhysicsSystem physicsSystem;

    physicsScene.Initialize();
    physicsSystem.Initialize(&physicsScene);

    std::cout << "Physics System Initialized." << std::endl;

    // --------------------------------------------------------
    // 3. 创建场景实体
    // --------------------------------------------------------

    // A. 地面 (Static)
    {
        Entity ground = registry.create();
        
        auto& t = registry.emplace<TransformComponent>(ground);
        t.setPosition(Vector3(0, -5, 0)); // 地面下移
        t.setScale(Vector3(1, 1, 1));
        
        auto& c = registry.emplace<ColliderComponent>(ground);
        c.setType(ColliderType::Box);
        c.setSize(Vector3(50.0f, 2.0f, 50.0f)); // 一个巨大的扁盒子

        auto& rb = registry.emplace<RigidBodyComponent>(ground);
        rb.setMass(0.0f); // 质量0 = 静态
        rb.setFriction(0.8f);
        rb.setRestitution(0.8f);//弹性系数
    }

    // B. 动态方块 (Dynamic)
    // 生成一个 3x3x3 的方块阵列
    for (int y = 0; y < 1; y++) {
        for (int x = 0; x < 3; x++) {
            Entity cube = registry.create();
            
            auto& t = registry.emplace<TransformComponent>(cube);
            // 随机一点位置，让它们掉落得自然些
            t.setPosition(Vector3(x * 1.5f - 1.0f, 10.0f + y * 2.0f, 0)); 
            t.setRotation(Quaternion::IDENTITY);
            t.setScale(Vector3(1, 1, 1));

            auto& c = registry.emplace<ColliderComponent>(cube);
            c.setType(ColliderType::Box);
            c.setSize(Vector3(1.0f,1.0f,1.0f));

            auto& rb = registry.emplace<RigidBodyComponent>(cube);
            rb.setMass(1.0f); // 质量1 = 动态
            rb.setFriction(0.1f);
            rb.setRestitution(0.5f);
        }
    }

    // --------------------------------------------------------
    // 4. 主循环
    // --------------------------------------------------------
    float lastTime = 0.0f;

    while (!glfwWindowShouldClose(window)) {
        float currentTime = (float)glfwGetTime();
        float deltaTime = currentTime - lastTime;
        lastTime = currentTime;

        // A. 物理 Tick
        physicsSystem.Tick(deltaTime, registry);

        // B. 渲染设置
        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // 设置简单的 透视投影 和 摄像机 (模拟 gluPerspective / gluLookAt)
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        // 简单的透视矩阵参数: fov=45, aspect, near=0.1, far=100
        float aspect = (float)WIN_WIDTH / (float)WIN_HEIGHT;
        float fH = tan(45.0f / 360.0f * 3.14159f) * 0.1f;
        float fW = fH * aspect;
        glFrustum(-fW, fW, -fH, fH, 0.1f, 100.0f);

        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        // 相机位置 (0, 10, 20) 看向 (0, 0, 0)
        glTranslatef(0.0f, -5.0f, -30.0f); 
        glRotatef(20.0f, 1.0f, 0.0f, 0.0f); // 稍微向下看

        // C. 遍历 ECS 并渲染
        auto view = registry.view<TransformComponent, ColliderComponent>();
        
        for (auto entity : view) {
            auto& t = view.get<TransformComponent>(entity);
            auto& c = view.get<ColliderComponent>(entity);

            glPushMatrix();
            
            // 1. 应用变换 (Pos + Rot + Scale)
            ApplyTransform(t);

            // 2. 加上 Collider 的 Offset (如果有)
            glTranslatef(c.getOffset().x, c.getOffset().y, c.getOffset().z);

            // 3. 根据类型画图
            if (c.getType() == ColliderType::Box) {
                // 地面画绿色，方块画橙色
                if (c.getSize().x > 10.0f) 
                    DrawCube(Vector3(0.2f, 0.8f, 0.2f), c.getSize());
                else 
                    DrawCube(Vector3(1.0f, 0.5f, 0.2f), c.getSize());
            }
            // (Sphere 和 Capsule 暂时用 Box 代替显示，因为画球比较麻烦)
            
            glPopMatrix();
        }

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    // 5. 清理
    physicsSystem.Shutdown();
    glfwTerminate();

    return 0;
}
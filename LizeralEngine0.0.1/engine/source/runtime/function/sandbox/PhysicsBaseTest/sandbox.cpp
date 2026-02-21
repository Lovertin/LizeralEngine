#include <iostream>
#include <vector>
#include <cmath>
#include <memory>

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

// 3. 射线与命中数据结构 (确保你已经创建了这两个文件)
#include "runtime/function/physics/phsicsEntityheaders.h"

// 4. 相机与输入系统头文件
#include "runtime/function/input/input.h"
#include "runtime/function/ecs/components/Camera/CameraComponent.h"
#include "runtime/function/ecs/components/Camera/CameraControlComponent.h"
#include "runtime/function/render/CameraControlSystem/CameraControlSystem.h"
#include "runtime/function/render/CameraSystem/CameraSystem.h"
#include "runtime/resource/resourceManager/resourceManager.h"
#include "runtime/function/pipeline/mesh/mesh.h"
#include "runtime/function/pipeline/shader/shader.h"
#include "runtime/function/pipeline/Material/Material.h"
#include "runtime/function/pipeline/Material/PBRMaterial.h"

using namespace Lizeral;

// --- 全局变量 ---
const int WIN_WIDTH = 1280;
const int WIN_HEIGHT = 720;

// --- 辅助：将 Lizeral Matrix4x4 转换为 OpenGL 格式并加载 ---
void LoadEngineMatrixToOpenGL(GLenum mode, const Matrix4x4& m) {
    float glMat[16];
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

void DrawMeshLegacy(const std::shared_ptr<Mesh>& mesh, const Vector3& color) {
    if (!mesh) return;

    const auto& verts = mesh->GetVertices();
    const auto& indices = mesh->GetIndices();

    glColor3f(color.x, color.y, color.z);
    
    // 禁用剔除，保证正反面都能看到
    glDisable(GL_CULL_FACE); 

    glBegin(GL_TRIANGLES);
    if (!indices.empty()) {
        // 如果有索引数据
        for (unsigned int idx : indices) {
            const auto& v = verts[idx];
            glNormal3f(v.Normal.x, v.Normal.y, v.Normal.z);
            glVertex3f(v.Position.x, v.Position.y, v.Position.z);
        }
    } else {
        // 如果没有索引（直接平铺的顶点）
        for (const auto& v : verts) {
            glNormal3f(v.Normal.x, v.Normal.y, v.Normal.z);
            glVertex3f(v.Position.x, v.Position.y, v.Position.z);
        }
    }
    glEnd();
    
    glEnable(GL_CULL_FACE); // 恢复状态
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

// --- 辅助：画一个简单的立方体 (新增边缘高亮) ---
void DrawCube(const Vector3& color, const Vector3& size, bool drawEdges = false) {
    float hx = size.x * 0.5f;
    float hy = size.y * 0.5f;
    float hz = size.z * 0.5f;

    // 1. 开启 Polygon Offset，防止线框和实心面深度冲突 (Z-fighting)
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1.0f, 1.0f);

    // 画实心面
    glColor3f(color.x, color.y, color.z);
    glBegin(GL_QUADS);
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

    glDisable(GL_POLYGON_OFFSET_FILL);

    // 2. 画黑色的边缘线
    if (drawEdges) {
        glColor3f(0.0f, 0.0f, 0.0f); // 黑色边框
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); // 切换为线框模式
        glLineWidth(2.0f); // 加粗一点线条

        glBegin(GL_QUADS);
        glVertex3f(-hx, -hy,  hz); glVertex3f( hx, -hy,  hz); glVertex3f( hx,  hy,  hz); glVertex3f(-hx,  hy,  hz);
        glVertex3f(-hx, -hy, -hz); glVertex3f(-hx,  hy, -hz); glVertex3f( hx,  hy, -hz); glVertex3f( hx, -hy, -hz);
        glVertex3f(-hx,  hy, -hz); glVertex3f(-hx,  hy,  hz); glVertex3f( hx,  hy,  hz); glVertex3f( hx,  hy, -hz);
        glVertex3f(-hx, -hy, -hz); glVertex3f( hx, -hy, -hz); glVertex3f( hx, -hy,  hz); glVertex3f(-hx, -hy,  hz);
        glVertex3f( hx, -hy, -hz); glVertex3f( hx,  hy, -hz); glVertex3f( hx,  hy,  hz); glVertex3f( hx, -hy,  hz);
        glVertex3f(-hx, -hy, -hz); glVertex3f(-hx, -hy,  hz); glVertex3f(-hx,  hy,  hz); glVertex3f(-hx,  hy, -hz);
        glEnd();

        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); // 恢复为填充模式
    }
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
    
    GLFWwindow* window = glfwCreateWindow(WIN_WIDTH, WIN_HEIGHT, "Lizeral Engine - Raycast Sandbox", NULL, NULL);
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
    
    PhysicsScene physicsScene;
    PhysicsSystem physicsSystem;
    physicsScene.Initialize();
    physicsSystem.Initialize(&physicsScene);

    Input::GetInstance().Init(window);
    
    // 【建议】：为了实现“右键按下才漫游”，请确保你在 Input::Init() 里
    // 暂时注释掉 `glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);`
    // 让初始状态是普通光标显示，否则你看不见鼠标点哪里。

    CameraSystem cameraSystem;
    CameraControlSystem cameraControlSystem;

    std::cout << "All Systems Initialized." << std::endl;
    std::cout << "[Controls] Hold Right Mouse Button + WASD/QE to Roam." << std::endl;
    std::cout << "[Controls] Click Left Mouse Button to select a cube." << std::endl;

    ResourceManager::GetInstance().SetRootPath("");
    auto testDolphinMesh = ResourceManager::GetInstance().Load<Mesh>("C:\\Lizeral Engine\\LizeralEngine0.0.1\\asset\\model\\dolphinLowPoly.model");

    TransformComponent dolphinTrans;
    dolphinTrans.setPosition(Vector3(0.0f, 15.0f, 0.0f)); // 悬浮在方块上方 15 米处
    dolphinTrans.setScale(Vector3(1.0f, 1.0f, 1.0f));

    auto pbrShader = std::make_shared<Shader>(
            "C:\\Lizeral Engine\\LizeralEngine0.0.1\\engine\\source\\runtime\\function\\sandbox\\shaderTest\\solid.vert", 
            "C:\\Lizeral Engine\\LizeralEngine0.0.1\\engine\\source\\runtime\\function\\sandbox\\shaderTest\\solid.frag"
        );

    auto dolphinPBR = std::make_shared<PBRMaterial>(pbrShader);

    dolphinPBR->m_Albedo = Vector3(1.0f, 0.86f, 0.57f); // 黄金的物理反照率
    dolphinPBR->m_Metallic = 1.0f;                      // 100% 纯金属
    dolphinPBR->m_Roughness = 0.2f;                     // 表面非常光滑，只有微小划痕
    dolphinPBR->m_AO = 1.0f;


    // --------------------------------------------------------
    // 3. 创建场景实体
    // --------------------------------------------------------

    // --- A. 创建主相机 ---
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

    // --- B. 地面 ---
    {
        Entity ground = registry.create();
        auto& t = registry.emplace<TransformComponent>(ground);
        t.setPosition(Vector3(0, -10, 0));
        t.setScale(Vector3(1, 1, 1));
        
        auto& c = registry.emplace<ColliderComponent>(ground);
        c.setType(ColliderType::Box);
        c.setSize(Vector3(100.0f, 2.0f, 100.0f));

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

        // --- 1. 输入处理 ---
        Input::GetInstance().Tick();
        if (Input::GetInstance().GetKey(Key::ESC)) glfwSetWindowShouldClose(window, true);

        // --- 2. 射线检测 (Raycast) 核心逻辑 ---
        if (Input::GetInstance().GetMouseButtonDown(MouseButton::Left)) {
            // 获取相机组件数据
            auto& mainCam = registry.get<CameraComponent>(cameraEntity);
            auto& camTrans = registry.get<TransformComponent>(cameraEntity);

            Vector2 mousePos = Input::GetInstance().GetMousePosition();
            
            // a) 将屏幕坐标转换为 NDC 坐标 [-1, 1]
            float ndcX = (2.0f * mousePos.x) / WIN_WIDTH - 1.0f;
            float ndcY = 1.0f - (2.0f * mousePos.y) / WIN_HEIGHT; // Y轴向下，需要翻转

            // b) 准备矩阵求逆 (假设你的数学库支持 Matrix4x4 inverse)
            Matrix4x4 invProj = mainCam.getProjectionMatrix().inverse();
            Matrix4x4 invView = mainCam.getViewMatrix().inverse();

            // c) 计算射线在 View 空间的坐标
            // 这里我们用一种无需 Vector4 类的方法，直接利用矩阵乘法解算远近两个点
            // 点1：近平面
            Vector3 nearPointNDC(ndcX, ndcY, -1.0f);
            Vector3 nearPointView = invProj * nearPointNDC; // 假设重载了 operator*(Vector3) 做齐次除法
            Vector3 nearPointWorld = invView * nearPointView;

            // 点2：远平面
            Vector3 farPointNDC(ndcX, ndcY, 1.0f);
            Vector3 farPointView = invProj * farPointNDC;
            Vector3 farPointWorld = invView * farPointView;

            // 计算方向
            Vector3 rayDir = farPointWorld - nearPointWorld;
            rayDir.normalise();

            // 发射射线！
            Vector3 rayOrigin = camTrans.getPosition(); // 相机位置作为起点
            Ray ray(rayOrigin, rayDir);
            RaycastHit hitInfo;

            // 调用物理引擎进行检测
            if (physicsSystem.Raycast(ray, hitInfo)) {
                std::cout << "\n===============================" << std::endl;
                std::cout << "[Raycast Hit!] Target Entity ID : " << hitInfo.entity << std::endl;
                std::cout << "   Hit Position : (" << hitInfo.point.x << ", " << hitInfo.point.y << ", " << hitInfo.point.z << ")" << std::endl;
                std::cout << "   Hit Normal   : (" << hitInfo.normal.x << ", " << hitInfo.normal.y << ", " << hitInfo.normal.z << ")" << std::endl;
                std::cout << "===============================\n" << std::endl;
            } else {
                std::cout << "[Raycast] Missed!" << std::endl;
            }
        }

        // --- 3. 物理与系统模拟 ---
        physicsSystem.Tick(deltaTime, registry);

        // 相机控制 (右键漫游)
        cameraControlSystem.Tick(registry, window); 
        // 重新计算相机 View/Proj 矩阵
        cameraSystem.Tick(registry);

        // --- 4. 渲染 ---
        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // 绑定相机矩阵
        auto& mainCam = registry.get<CameraComponent>(cameraEntity);
        auto& camTrans = registry.get<TransformComponent>(cameraEntity);

        LoadEngineMatrixToOpenGL(GL_PROJECTION, mainCam.getProjectionMatrix());
        LoadEngineMatrixToOpenGL(GL_MODELVIEW, mainCam.getViewMatrix());

        // 绘制物理实体
        auto view = registry.view<TransformComponent, ColliderComponent>();
        
        for (auto entity : view) {
            auto& t = view.get<TransformComponent>(entity);
            auto& c = view.get<ColliderComponent>(entity);

            glPushMatrix();
            
            ApplyTransform(t);
            glTranslatef(c.getOffset().x, c.getOffset().y, c.getOffset().z);

            // 绘制区分：大物体(地面)不画边框，小方块画黑色线框
            if (c.getType() == ColliderType::Box) {
                if (c.getSize().x > 10.0f) {
                    DrawCube(Vector3(0.2f, 0.8f, 0.2f), c.getSize(), false);
                } else {
                    // 动态小方块画上黑边，更方便视觉定位
                    DrawCube(Vector3(1.0f, 0.5f, 0.2f), c.getSize(), true);
                }
            }
            
            glPopMatrix();
        }

        if (testDolphinMesh && pbrShader) {
            // 让你在沙盒里也能用上下键动态缩放海豚
            if (Input::GetInstance().GetKey(Key::UP)) dolphinTrans.setScale(dolphinTrans.getScale() * 1.05f);
            if (Input::GetInstance().GetKey(Key::DOWN)) dolphinTrans.setScale(dolphinTrans.getScale() * 0.95f);

            // 让海豚自转
            // static float rotationAngle = 0.0f;
            // rotationAngle += deltaTime * 50.0f; 
            // dolphinTrans.setRotation(Quaternion(Vector3(0, 1, 0), rotationAngle * (3.14159f / 180.0f)));

            // ==========================================
            // 现代管线接管开始！
            // ==========================================
            pbrShader->Bind();

            // 1. 系统级矩阵与相机数据
            pbrShader->SetUniformMat4f("model", dolphinTrans.getMatrix());
            pbrShader->SetUniformMat4f("view", mainCam.getViewMatrix());
            pbrShader->SetUniformMat4f("projection", mainCam.getProjectionMatrix());
            
            // 【关键】：PBR 极度依赖观察视角，必须把相机的世界坐标传进去！
            pbrShader->SetUniformVector3f("u_camPos", camTrans.getPosition());

            // 2. 光源数据
            Vector3 lightDir(1.0f, 1.0f, 1.0f); 
            lightDir.normalise();
            pbrShader->SetUniformVector3f("u_lightDir", lightDir);
            pbrShader->SetUniformVector3f("u_lightColor", Vector3(3.0f, 3.0f, 3.0f)); // 稍微加强光强，体现 HDR

            // 3. 材质数据：直接调用接口！引擎不管它是 PBR 还是普通的
            dolphinPBR->BindAndApply();

            glDisable(GL_CULL_FACE);
            testDolphinMesh->Draw();
            glEnable(GL_CULL_FACE);

            glUseProgram(0);
        }

        glfwSwapBuffers(window);
    }

    // 5. 清理
    physicsSystem.Shutdown();
    glfwTerminate();

    return 0;
}
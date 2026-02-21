#include <iostream>
#include <memory>
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <cmath>

#include "runtime/function/ecs/registry.h"
#include "runtime/function/ecs/components/Transform/TransformComponent.h"
#include "runtime/function/ecs/components/Camera/CameraComponent.h"
#include "runtime/resource/resourceManager/resourceManager.h"
#include "runtime/function/res_type/mesh/mesh.h"
#include "runtime/function/res_type/shader/shader.h"

using namespace Lizeral;

// =======================================================
// 全局变量：漫游相机状态
// =======================================================
Vector3 g_cameraPos(0.0f, 2.0f, 10.0f);   // 相机初始位置
Vector3 g_cameraFront(0.0f, 0.0f, -1.0f); // 相机正前方方向
Vector3 g_cameraUp(0.0f, 1.0f, 0.0f);     // 世界的正上方

float g_deltaTime = 0.0f; // 当前帧与上一帧的时间差
float g_lastFrame = 0.0f; // 上一帧的时间

// 鼠标控制变量
bool g_firstMouse = true;
float g_yaw   = -90.0f; // 初始偏航角（指向 -Z 轴）
float g_pitch =  0.0f;
float g_lastX =  1280.0f / 2.0;
float g_lastY =  720.0 / 2.0;

// =======================================================
// 辅助数学函数：构建 LookAt 视图矩阵
// =======================================================
Matrix4x4 CreateLookAtMatrix(const Vector3& eye, const Vector3& target, const Vector3& up) {
    // 1. 计算前向向量 (Forward)
    Vector3 f = target - eye;
    f.normalise();
    
    // 2. 计算右向量 (Right) = Forward x Up
    Vector3 s(f.y * up.z - f.z * up.y, f.z * up.x - f.x * up.z, f.x * up.y - f.y * up.x);
    s.normalise();
    
    // 3. 计算真正的上向量 (Up) = Right x Forward
    Vector3 u(s.y * f.z - s.z * f.y, s.z * f.x - s.x * f.z, s.x * f.y - s.y * f.x);

    // 4. 构建 4x4 矩阵
    Matrix4x4 mat;
    // 第一行 (Right)
    mat.m_mat[0][0] = s.x; mat.m_mat[0][1] = s.y; mat.m_mat[0][2] = s.z; 
    mat.m_mat[0][3] = -(s.x * eye.x + s.y * eye.y + s.z * eye.z);
    // 第二行 (Up)
    mat.m_mat[1][0] = u.x; mat.m_mat[1][1] = u.y; mat.m_mat[1][2] = u.z; 
    mat.m_mat[1][3] = -(u.x * eye.x + u.y * eye.y + u.z * eye.z);
    // 第三行 (-Forward)
    mat.m_mat[2][0] = -f.x; mat.m_mat[2][1] = -f.y; mat.m_mat[2][2] = -f.z; 
    mat.m_mat[2][3] = (f.x * eye.x + f.y * eye.y + f.z * eye.z);
    // 第四行
    mat.m_mat[3][0] = 0.0f; mat.m_mat[3][1] = 0.0f; mat.m_mat[3][2] = 0.0f; 
    mat.m_mat[3][3] = 1.0f;

    return mat;
}

// =======================================================
// 回调与输入处理
// =======================================================
void mouse_callback(GLFWwindow* window, double xposIn, double yposIn) {
    float xpos = static_cast<float>(xposIn);
    float ypos = static_cast<float>(yposIn);

    if (g_firstMouse) {
        g_lastX = xpos;
        g_lastY = ypos;
        g_firstMouse = false;
    }

    float xoffset = xpos - g_lastX;
    float yoffset = g_lastY - ypos; // 注意这里是反过来的
    g_lastX = xpos;
    g_lastY = ypos;

    float sensitivity = 0.1f;
    xoffset *= sensitivity;
    yoffset *= sensitivity;

    g_yaw   += xoffset;
    g_pitch += yoffset;

    if (g_pitch > 89.0f)  g_pitch = 89.0f;
    if (g_pitch < -89.0f) g_pitch = -89.0f;

    Vector3 front;
    front.x = cos(g_yaw * (3.14159f / 180.0f)) * cos(g_pitch * (3.14159f / 180.0f));
    front.y = sin(g_pitch * (3.14159f / 180.0f));
    front.z = sin(g_yaw * (3.14159f / 180.0f)) * cos(g_pitch * (3.14159f / 180.0f));
    front.normalise();
    g_cameraFront = front;
}

// 处理 WASD 移动和海豚缩放
void processInput(GLFWwindow* window, TransformComponent& dolphinTrans) {
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);

    float cameraSpeed = 5.0f * g_deltaTime;
    if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS)
        cameraSpeed *= 3.0f; // 按住 Shift 加速飞行

    // 相机移动 (WASD + QE)
    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
        g_cameraPos = g_cameraPos + (g_cameraFront * cameraSpeed);
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
        g_cameraPos = g_cameraPos - (g_cameraFront * cameraSpeed);
    
    // 计算右向量并移动
    Vector3 right(g_cameraFront.y * g_cameraUp.z - g_cameraFront.z * g_cameraUp.y, 
                  g_cameraFront.z * g_cameraUp.x - g_cameraFront.x * g_cameraUp.z, 
                  g_cameraFront.x * g_cameraUp.y - g_cameraFront.y * g_cameraUp.x);
    right.normalise();
    
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
        g_cameraPos = g_cameraPos - (right * cameraSpeed);
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
        g_cameraPos = g_cameraPos + (right * cameraSpeed);
    if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS)
        g_cameraPos = g_cameraPos + (g_cameraUp * cameraSpeed);
    if (glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS)
        g_cameraPos = g_cameraPos - (g_cameraUp * cameraSpeed);

    // ==========================================
    // 关键排查功能：方向键实时调整海豚大小！
    // ==========================================
    float scaleSpeed = 1.05f; // 每次放大 5%
    if (glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS) {
        Vector3 curScale = dolphinTrans.getScale();
        dolphinTrans.setScale(curScale * scaleSpeed);
    }
    if (glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS) {
        Vector3 curScale = dolphinTrans.getScale();
        dolphinTrans.setScale(curScale * (1.0f / scaleSpeed));
    }
}

// =======================================================
// 材质与组件定义
// =======================================================
class SolidMaterial {
public:
    std::shared_ptr<Shader> shader;
    Vector3 albedo {0.2f, 0.6f, 0.9f}; // 海豚蓝
    Vector3 lightDir {1.0f, 1.0f, 1.0f}; // 阳光从右上方打来
    Vector3 lightColor {1.0f, 1.0f, 1.0f}; // 白光

    SolidMaterial(std::shared_ptr<Shader> s) : shader(s) {}

    void BindAndApply() {
        if (!shader) return;
        shader->Bind();
        shader->SetUniformVector3f("u_albedo", albedo);
        shader->SetUniformVector3f("u_lightDir", lightDir);
        shader->SetUniformVector3f("u_lightColor", lightColor);
    }
};

class MeshRendererComponent : public Component {
public:
    std::shared_ptr<Mesh> mesh;
    std::shared_ptr<SolidMaterial> material;
    MeshRendererComponent(std::shared_ptr<Mesh> m, std::shared_ptr<SolidMaterial> mat) 
        : mesh(m), material(mat) {}
};

// =======================================================
// 主循环
// =======================================================
int main() {
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    GLFWwindow* window = glfwCreateWindow(1280, 720, "Lizeral Engine - Dolphin Renderer", NULL, NULL);
    glfwMakeContextCurrent(window);
    gladLoadGLLoader((GLADloadproc)glfwGetProcAddress);

    // 捕获鼠标，隐藏光标用于沉浸式漫游
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
    glfwSetCursorPosCallback(window, mouse_callback);

    // 开启深度测试
    glEnable(GL_DEPTH_TEST);
    // 【关键排查】：禁用背面剔除，防止模型法线装反导致隐形
    glDisable(GL_CULL_FACE);

    ResourceManager::GetInstance().SetRootPath("");
    

    std::cout << "Loading resources..." << std::endl;
    // 替换为你电脑上的绝对路径
    auto dolphinMesh = ResourceManager::GetInstance().Load<Mesh>("C:\\Lizeral Engine\\LizeralEngine0.0.1\\asset\\model\\dolphinLowPoly.model");
    
    auto solidShader = std::make_shared<Shader>(
        "C:\\Lizeral Engine\\LizeralEngine0.0.1\\engine\\source\\runtime\\function\\sandbox\\shaderTest\\solid.vert", 
        "C:\\Lizeral Engine\\LizeralEngine0.0.1\\engine\\source\\runtime\\function\\sandbox\\shaderTest\\solid.frag"
    );

    auto dolphinMaterial = std::make_shared<SolidMaterial>(solidShader);

    Registry registry;

    // 创建相机组件
    Entity camera = registry.create();
    auto& camTrans = registry.emplace<TransformComponent>(camera);
    auto& camComp = registry.emplace<CameraComponent>(camera);
    camComp.setFov(45.0f);
    camComp.setAspect(1280.0f / 720.0f);
    camComp.setzNear(0.1f);
    camComp.setzFar(1000.0f); // 把远裁剪面调大，防止模型太大被切掉
    // 初始投影矩阵
    camComp.setProjectionMatrix(camComp.buildPerspective(45.0f * (3.14159f / 180.0f), 1280.0f/720.0f, 0.1f, 1000.0f));

    // 创建海豚实体
    Entity dolphin = registry.create();
    auto& dolTrans = registry.emplace<TransformComponent>(dolphin);
    dolTrans.setPosition(Vector3(0.0f, 0.0f, 0.0f));
    dolTrans.setScale(Vector3(1.0f, 1.0f, 1.0f));
    registry.emplace<MeshRendererComponent>(dolphin, dolphinMesh, dolphinMaterial);

    std::cout << "Ready to Render! Use WASD to fly, Mouse to look, UP/DOWN arrows to scale the dolphin." << std::endl;

    while (!glfwWindowShouldClose(window)) {
        float currentFrame = glfwGetTime();
        g_deltaTime = currentFrame - g_lastFrame;
        g_lastFrame = currentFrame;

        // 1. 处理输入（包含对海豚的缩放）
        processInput(window, dolTrans);

        // 2. 更新相机视图矩阵 (核心！)
        camTrans.setPosition(g_cameraPos);
        Matrix4x4 viewMat = CreateLookAtMatrix(g_cameraPos, g_cameraPos + g_cameraFront, g_cameraUp);
        camComp.setViewMatrix(viewMat);
        Matrix4x4 projMat = camComp.getProjectionMatrix();

        // 3. 渲染指令
        glClearColor(0.15f, 0.15f, 0.15f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        auto view = registry.view<TransformComponent, MeshRendererComponent>();
        for (auto entity : view) {
            auto& trans = view.get<TransformComponent>(entity);
            auto& renderer = view.get<MeshRendererComponent>(entity);

            if (renderer.mesh && renderer.material) {
                renderer.material->BindAndApply();
                
                Matrix4x4 modelMat = trans.getMatrix(); // 获取模型矩阵
                renderer.material->shader->SetUniformMat4f("model", modelMat);
                renderer.material->shader->SetUniformMat4f("view", viewMat);
                renderer.material->shader->SetUniformMat4f("projection", projMat);

                renderer.mesh->Draw();
            }
        }

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwTerminate();
    return 0;
}
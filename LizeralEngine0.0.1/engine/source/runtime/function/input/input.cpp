#include "input.h"
#include <cstring> // for memcpy, memset

namespace Lizeral {

    // --- 构造函数 (初始化成员变量) ---
    Input::Input() {
        std::memset(m_keys, 0, sizeof(m_keys));
        std::memset(m_lastKeys, 0, sizeof(m_lastKeys));
        std::memset(m_mouseButtons, 0, sizeof(m_mouseButtons));
        std::memset(m_lastMouseButtons, 0, sizeof(m_lastMouseButtons));
    }

    void Input::Init(GLFWwindow* window) {
        // 绑定回调 (这些回调函数仍然是 static 的，这是必须的)
        glfwSetKeyCallback(window, KeyCallback);
        glfwSetMouseButtonCallback(window, MouseButtonCallback);
        glfwSetScrollCallback(window, ScrollCallback);
        glfwSetCursorPosCallback(window, CursorPosCallback);
        glfwSetWindowFocusCallback(window, WindowFocusCallback);
        
        // 禁用系统光标，开启 FPS 模式
        glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
    }

    void Input::Tick() {
        // 1. 在 PollEvents 之前，清空上一帧的瞬时状态 (Delta)
        // 必须在这里清零，因为接下来的 glfwPollEvents 会调用 CursorPosCallback 进行累加
        m_mouseDelta = { 0, 0 };
        m_scrollDelta = 0.0f;

        // 2. 将当前帧状态复制到上一帧 (为下一帧的 GetKeyDown/Up 做准备)
        std::memcpy(m_lastKeys, m_keys, sizeof(m_keys));
        std::memcpy(m_lastMouseButtons, m_mouseButtons, sizeof(m_mouseButtons));

        // 3. 让 GLFW 处理系统消息，这将触发 Callbacks 并更新 m_keys / m_mouseDelta
        glfwPollEvents();
    }

    void Input::ResetMouse() {
        s_firstMouseInput = true;
        m_mouseDelta = Vector2(0, 0); 
    }

    // --- 静态回调实现 (Static Callbacks) ---
    // 必须通过 GetInstance() 访问单例的成员变量

    void Input::KeyCallback(GLFWwindow* window, int key, int scancode, int action, int mods) {
        if (key >= 0 && key < GLFW_KEY_LAST) {
            GetInstance().m_keys[key] = (action != GLFW_RELEASE);
        }
    }

    void Input::CursorPosCallback(GLFWwindow* window, double xpos, double ypos) {
        Input& self = GetInstance(); // 获取实例引用

        Vector2 newPos((float)xpos, (float)ypos);
        
        if (self.s_firstMouseInput) {
            self.m_mousePos = newPos;
            self.s_firstMouseInput = false;
            return;
        }

        // 累加 Delta (防止一帧内多次回调覆盖)
        self.m_mouseDelta += newPos - self.m_mousePos;
        self.m_mousePos = newPos;
    }

    void Input::MouseButtonCallback(GLFWwindow* window, int button, int action, int mods) {
        if (button >= 0 && button < GLFW_MOUSE_BUTTON_LAST) {
            GetInstance().m_mouseButtons[button] = (action != GLFW_RELEASE);
        }
    }

    void Input::ScrollCallback(GLFWwindow* window, double xoffset, double yoffset) {
        GetInstance().m_scrollDelta = (float)yoffset;
    }

    void Input::WindowFocusCallback(GLFWwindow* window, int focused) {
        if (focused) {
            // 窗口获得焦点：重置鼠标状态，重新锁定
            GetInstance().ResetMouse();
            glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
        } else {
            // 窗口失去焦点：显示光标
            glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
        }
    }

    // --- 成员函数实现 (Member Functions) ---
    // 注意：这里去掉了 static 关键字

    bool Input::GetKey(Key key) {
        return m_keys[(int)key];
    }

    bool Input::GetKeyDown(Key key) {
        return m_keys[(int)key] && !m_lastKeys[(int)key];
    }

    bool Input::GetKeyUp(Key key) {
        return !m_keys[(int)key] && m_lastKeys[(int)key];
    }

    bool Input::GetMouseButton(MouseButton button) {
        return m_mouseButtons[(int)button];
    }

    bool Input::GetMouseButtonDown(MouseButton button) {
        return m_mouseButtons[(int)button] && !m_lastMouseButtons[(int)button];
    }

    bool Input::GetMouseButtonUp(MouseButton button) {
        return !m_mouseButtons[(int)button] && m_lastMouseButtons[(int)button];
    }

    Vector2 Input::GetMousePosition() {
        return m_mousePos;
    }

    float Input::GetMouseScroll() {
        return m_scrollDelta;
    }

    Vector2 Input::GetMouseDelta() {
        return m_mouseDelta;
    }

} // namespace Lizeral
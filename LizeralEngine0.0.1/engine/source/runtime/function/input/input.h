#pragma once
#include "runtime/core/math/vector2.h"
#include <GLFW/glfw3.h>

namespace Lizeral{
    enum class Key{
        W=GLFW_KEY_W,
        A=GLFW_KEY_A,
        S=GLFW_KEY_S,
        D=GLFW_KEY_D,
        E=GLFW_KEY_E,
        Q=GLFW_KEY_Q,
        ESC=GLFW_KEY_ESCAPE,
        SPACE=GLFW_KEY_SPACE,
        R=GLFW_KEY_R,
        LEFT_SHIFT=GLFW_KEY_LEFT_SHIFT
    };

    enum class MouseButton {
        Left = GLFW_MOUSE_BUTTON_LEFT,
        Right = GLFW_MOUSE_BUTTON_RIGHT,
        Middle = GLFW_MOUSE_BUTTON_MIDDLE
    };

    class Input{
    public:
        // 删除拷贝构造和赋值操作
        Input(const Input&) = delete;
        Input& operator=(const Input&) = delete;

        // 获取单例实例
        static Input& GetInstance() {
            static Input instance;
            return instance;
        }

        void Init(GLFWwindow* window);
        void Tick();

        // 键盘
        bool GetKey(Key key);
        bool GetKeyDown(Key key);
        bool GetKeyUp(Key key);

        // 鼠标
        bool GetMouseButton(MouseButton button);
        bool GetMouseButtonDown(MouseButton button);
        bool GetMouseButtonUp(MouseButton button);
        
        Vector2 GetMousePosition();
        float GetMouseScroll();

        // 核心：这一帧鼠标移动了多少 (dx, dy)
        Vector2 GetMouseDelta();

        //重置鼠标状态 
        void ResetMouse();


    private:
        //私有构造函数
        Input();
        ~Input() = default;

        // GLFW 回调函数 (必须是 static)
        static void KeyCallback(GLFWwindow* window, int key, int scancode, int action, int mods);
        static void MouseButtonCallback(GLFWwindow* window, int button, int action, int mods);
        static void ScrollCallback(GLFWwindow* window, double xoffset, double yoffset);
        static void CursorPosCallback(GLFWwindow* window, double xpos, double ypos);
        //回调查询是否选中当前窗口
        static void WindowFocusCallback(GLFWwindow* window, int focused);

    private:
        bool m_keys[GLFW_KEY_LAST];      // 当前帧键盘
        bool m_lastKeys[GLFW_KEY_LAST];  // 上一帧键盘
        
        bool m_mouseButtons[GLFW_MOUSE_BUTTON_LAST];
        bool m_lastMouseButtons[GLFW_MOUSE_BUTTON_LAST];

        Vector2 m_mousePos;
        Vector2 m_mouseDelta;
        float m_scrollDelta;

        bool s_firstMouseInput;
    };
}

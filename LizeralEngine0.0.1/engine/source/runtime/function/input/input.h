#pragma once
#include "runtime/core/math/vector2.h"
#include <GLFW/glfw3.h>

namespace Lizeral{
    enum class Key{
        W=GLFW_KEY_W,
        A=GLFW_KEY_A,
        S=GLFW_KEY_S,
        D=GLFW_KEY_D,
        ESC=GLFW_KEY_ESCAPE,
        SPACE=GLFW_KEY_SPACE,
        R=GLFW_KEY_R
    };

    enum class MouseButton {
        Left = GLFW_MOUSE_BUTTON_LEFT,
        Right = GLFW_MOUSE_BUTTON_RIGHT,
        Middle = GLFW_MOUSE_BUTTON_MIDDLE
    };

    class Input{
    public:
        void Init(GLFWwindow* window);
        static void Tick();

        // 键盘
        static bool GetKey(Key key);
        static bool GetKeyDown(Key key);
        static bool GetKeyUp(Key key);

        // 鼠标
        static bool GetMouseButton(MouseButton button);
        static bool GetMouseButtonDown(MouseButton button);
        static bool GetMouseButtonUp(MouseButton button);
        
        static Vector2 GetMousePosition();
        static float GetMouseScroll();

        // 核心：这一帧鼠标移动了多少 (dx, dy)
        static Vector2 GetMouseDelta();

        //重置鼠标状态 
        static void ResetMouse();


    private:
        // GLFW 回调函数 (必须是 static)
        static void KeyCallback(GLFWwindow* window, int key, int scancode, int action, int mods);
        static void MouseButtonCallback(GLFWwindow* window, int button, int action, int mods);
        static void ScrollCallback(GLFWwindow* window, double xoffset, double yoffset);
        static void CursorPosCallback(GLFWwindow* window, double xpos, double ypos);
        //回调查询是否选中当前窗口
        static void WindowFocusCallback(GLFWwindow* window, int focused);

    private:
        static bool m_keys[GLFW_KEY_LAST];      // 当前帧键盘
        static bool m_lastKeys[GLFW_KEY_LAST];  // 上一帧键盘
        
        static bool m_mouseButtons[GLFW_MOUSE_BUTTON_LAST];
        static bool m_lastMouseButtons[GLFW_MOUSE_BUTTON_LAST];

        static Vector2 m_mousePos;
        static Vector2 m_mouseDelta;
        static float m_scrollDelta;

        static bool s_firstMouseInput;
    };
}

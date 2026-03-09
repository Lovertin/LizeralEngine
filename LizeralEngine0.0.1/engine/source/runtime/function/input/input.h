#pragma once
#include "runtime/core/math/vector2.h"
// 【彻底删除】：#include <GLFW/glfw3.h>

namespace Lizeral{

    // 【新增】：定义我们自己引擎的最大按键和鼠标键数量
    constexpr int MAX_KEYS = 512;
    constexpr int MAX_MOUSE_BUTTONS = 8;

    // 【重构】：不再依赖 GLFW 宏，直接使用干净的枚举值
    enum class Key{
        W = 0, A, S, D, E, Q, 
        ESC, SPACE, R, LEFT_SHIFT, DOWN, UP
        // 未来如果需要更多按键，直接往下加即可
    };

    enum class MouseButton {
        Left = 0, 
        Right, 
        Middle
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

        // 【删除】：void Init(GLFWwindow* window);
        
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

        // 重置鼠标状态 
        void ResetMouse();

        // 被动接收模式的接口（专供 Qt 等外部窗口系统调用）
        void SetKeyDown(Key key, bool isDown);
        void SetMouseButtonDown(MouseButton button, bool isDown);
        void SetMousePosition(float x, float y);
        void SetMouseScroll(float delta);

    private:
        // 私有构造函数
        Input();
        ~Input() = default;

        // 【彻底删除】：所有静态的 GLFW Callbacks
        
    private:
        // 替换为我们自己的常量
        bool m_keys[MAX_KEYS] = {false};       
        bool m_lastKeys[MAX_KEYS] = {false};   
        
        bool m_mouseButtons[MAX_MOUSE_BUTTONS] = {false};
        bool m_lastMouseButtons[MAX_MOUSE_BUTTONS] = {false};

        Vector2 m_mousePos{0.0f, 0.0f};
        Vector2 m_mouseDelta{0.0f, 0.0f};
        float m_scrollDelta{0.0f};

        bool s_firstMouseInput{true};
    };
}
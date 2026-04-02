#pragma once
#include "runtime/core/math/vector2.h"
// 【彻底删除】：#include <GLFW/glfw3.h>

namespace Lizeral{

    constexpr int MAX_KEYS = 512;
    constexpr int MAX_MOUSE_BUTTONS = 8;

    enum class Key{
        W = 0, A, S, D, E, Q, 
        ESC, SPACE, R, LEFT_SHIFT, DOWN, UP
    };

    enum class MouseButton {
        Left = 0, 
        Right, 
        Middle
    };

    class Input{
    public:
        Input(const Input&) = delete;
        Input& operator=(const Input&) = delete;

        static Input& GetInstance() {
            static Input instance;
            return instance;
        }

        
        void Tick();

        bool GetKey(Key key);
        bool GetKeyDown(Key key);
        bool GetKeyUp(Key key);

        bool GetMouseButton(MouseButton button);
        bool GetMouseButtonDown(MouseButton button);
        bool GetMouseButtonUp(MouseButton button);
        
        Vector2 GetMousePosition();
        float GetMouseScroll();

        Vector2 GetMouseDelta();

        void ResetMouse();
        void ResetState();

        void SetKeyDown(Key key, bool isDown);
        void SetMouseButtonDown(MouseButton button, bool isDown);
        void SetMousePosition(float x, float y);
        void SetMouseScroll(float delta);

    private:
        Input();
        ~Input() = default;

        
    private:
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

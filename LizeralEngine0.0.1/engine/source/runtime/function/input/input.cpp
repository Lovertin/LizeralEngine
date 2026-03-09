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

    // 【彻底删除】：Init 函数及其所有的 glfwSet...Callback

    void Input::Tick() {
        // 1. 将当前帧状态复制到上一帧 (为下一帧的 GetKeyDown/Up 做准备)
        std::memcpy(m_lastKeys, m_keys, sizeof(m_keys));
        std::memcpy(m_lastMouseButtons, m_mouseButtons, sizeof(m_mouseButtons));

        // 2. 清空这一帧的增量数据 (Delta)
        // 注意：因为我们是依靠外部框架主动调用 Set 喂数据的（Push 模型），
        // 当帧逻辑跑完后，清空 Delta 准备迎接 Qt 的下一次输入。
        m_mouseDelta = Vector2(0.0f, 0.0f);
        m_scrollDelta = 0.0f;
    }

    void Input::ResetMouse() {
        s_firstMouseInput = true;
    }

    // 【彻底删除】：所有的 KeyCallback, CursorPosCallback, MouseButtonCallback 等静态回调函数

    // --- 成员函数实现 (Member Functions) ---

    bool Input::GetKey(Key key) {
        return m_keys[static_cast<int>(key)];
    }

    bool Input::GetKeyDown(Key key) {
        int index = static_cast<int>(key);
        return m_keys[index] && !m_lastKeys[index];
    }

    bool Input::GetKeyUp(Key key) {
        int index = static_cast<int>(key);
        return !m_keys[index] && m_lastKeys[index];
    }

    bool Input::GetMouseButton(MouseButton button) {
        return m_mouseButtons[static_cast<int>(button)];
    }

    bool Input::GetMouseButtonDown(MouseButton button) {
        int index = static_cast<int>(button);
        return m_mouseButtons[index] && !m_lastMouseButtons[index];
    }

    bool Input::GetMouseButtonUp(MouseButton button) {
        int index = static_cast<int>(button);
        return !m_mouseButtons[index] && m_lastMouseButtons[index];
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

    void Input::SetKeyDown(Key key, bool isDown) {
        int keyIndex = static_cast<int>(key);
        if (keyIndex >= 0 && keyIndex < MAX_KEYS) {
            m_keys[keyIndex] = isDown;
        }
    }

    void Input::SetMouseButtonDown(MouseButton button, bool isDown) {
        int btnIndex = static_cast<int>(button);
        if (btnIndex >= 0 && btnIndex < MAX_MOUSE_BUTTONS) {
            m_mouseButtons[btnIndex] = isDown;
        }
    }

    void Input::SetMousePosition(float x, float y) {
        if (s_firstMouseInput) {
            m_mousePos = Vector2(x, y);
            s_firstMouseInput = false;
        }
        // 计算 Delta（当前帧移动量）
        m_mouseDelta.x += x - m_mousePos.x;
        m_mouseDelta.y += y - m_mousePos.y; 
        
        // 更新当前位置
        m_mousePos = Vector2(x, y);
    }

    void Input::SetMouseScroll(float delta) {
        m_scrollDelta += delta;
    }

} // namespace Lizeral
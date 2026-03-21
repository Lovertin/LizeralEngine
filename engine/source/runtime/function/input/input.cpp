#include "input.h"
#include <cstring> // for memcpy, memset

namespace Lizeral {

    Input::Input() {
        std::memset(m_keys, 0, sizeof(m_keys));
        std::memset(m_lastKeys, 0, sizeof(m_lastKeys));
        std::memset(m_mouseButtons, 0, sizeof(m_mouseButtons));
        std::memset(m_lastMouseButtons, 0, sizeof(m_lastMouseButtons));
    }

    void Input::Tick() {
        // prepare GetKeyDown/Up 
        std::memcpy(m_lastKeys, m_keys, sizeof(m_keys));
        std::memcpy(m_lastMouseButtons, m_mouseButtons, sizeof(m_mouseButtons));

        m_mouseDelta = Vector2(0.0f, 0.0f);
        m_scrollDelta = 0.0f;
    }

    void Input::ResetMouse() {
        s_firstMouseInput = true;
    }

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
        // caculate delta
        m_mouseDelta.x += x - m_mousePos.x;
        m_mouseDelta.y += y - m_mousePos.y; 

        m_mousePos = Vector2(x, y);
    }

    void Input::SetMouseScroll(float delta) {
        m_scrollDelta += delta;
    }

} // namespace Lizeral
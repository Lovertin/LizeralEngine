#include "input.h"
#include <cstring>

namespace Lizeral{
    bool Input::m_keys[GLFW_KEY_LAST] = { 0 };
    bool Input::m_lastKeys[GLFW_KEY_LAST] = { 0 };
    bool Input::m_mouseButtons[GLFW_MOUSE_BUTTON_LAST] = { 0 };
    bool Input::m_lastMouseButtons[GLFW_MOUSE_BUTTON_LAST] = { 0 };
    Vector2 Input::m_mousePos = { 0, 0 };
    Vector2 Input::m_mouseDelta = { 0, 0 };
    float Input::m_scrollDelta = 0.0f;
    static bool s_firstMouseInput = true;

    void Input::Init(GLFWwindow* window) {
        // 绑定回调
        glfwSetKeyCallback(window, KeyCallback);
        glfwSetMouseButtonCallback(window, MouseButtonCallback);
        glfwSetScrollCallback(window, ScrollCallback);
        glfwSetCursorPosCallback(window, CursorPosCallback);

        glfwSetWindowFocusCallback(window, WindowFocusCallback);
        
        // 这一步很重要：禁用系统光标，实现 FPS 视角的无限制移动
        // 如果是编辑器模式，可能需要按住右键才隐藏，这里先做个简单的
        glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
    }

    void Input::Tick(){
        // 1. 将当前帧状态复制到上一帧 (为下一帧的 GetKeyDown/Up 做准备)
        // 注意：这里要在 Polling 之前做，还是之后做，取决于你的 Loop 顺序
        // 最佳实践：
        // Frame Start -> Input::Tick() -> Copy Current to Last -> Clear Deltas -> PollEvents
        
        memcpy(m_lastKeys, m_keys, sizeof(m_keys));
        memcpy(m_lastMouseButtons, m_mouseButtons, sizeof(m_mouseButtons));
        
        // 2. 清空 Delta 类数据 (因为 Delta 是瞬时的)
        m_mouseDelta = { 0, 0 }; // 修正：Delta 不能在这里清零，要靠 CursorPosCallback 计算
        // 其实 Delta 的计算比较特殊，通常是在 Tick 结尾清零，或者在 Callback 里累加
        // 最简单的做法：
        // Delta = CurrentPos - LastFramePos
        // 所以我们需要记录 m_lastMousePos
        
        m_scrollDelta = 0.0f; // 滚轮是瞬发事件，每帧重置

        // 3. 让 GLFW 处理系统消息，这将触发上面的 Callbacks 并更新 m_keys
        glfwPollEvents();
    }

    // --- 回调实现 ---

    void Input::KeyCallback(GLFWwindow* window, int key, int scancode, int action, int mods) {
        if (key >= 0 && key < GLFW_KEY_LAST) {
            m_keys[key] = (action != GLFW_RELEASE);
        }
    }

    void Input::CursorPosCallback(GLFWwindow* window, double xpos, double ypos) {
         Vector2 newPos((float)xpos, (float)ypos);
        if(s_firstMouseInput){
            m_mousePos=newPos;
            s_firstMouseInput=false;
            return;
        }
        // 计算 Delta
        
        // 【修正】使用 += 累加，防止一帧内多次回调导致数据覆盖
        Vector2 frameDelta = newPos - m_mousePos;
        m_mouseDelta += frameDelta; 
        
        m_mousePos = newPos;
    }


    void Input::MouseButtonCallback(GLFWwindow* window, int button, int action, int mods) {
        // 过滤掉我们不关心的按键 ID（防止数组越界）
        if (button >= 0 && button < GLFW_MOUSE_BUTTON_LAST) {
            // action == GLFW_PRESS   -> 存为 true
            // action == GLFW_RELEASE -> 存为 false
            m_mouseButtons[button] = (action != GLFW_RELEASE);
        }
    }


    void Input::ScrollCallback(GLFWwindow* window, double xoffset, double yoffset) {
        // 我们主要关心 yoffset (垂直滚动：推+1，拉-1)
        // 直接累加或者赋值给 m_scrollDelta
        m_scrollDelta = (float)yoffset;
    }
    
    void Input::WindowFocusCallback(GLFWwindow* window, int focused){
        if (focused) {
            // 窗口获得了焦点 (focused == GLFW_TRUE)
            // 意味着玩家切回来了，我们需要重置 "第一次输入" 的标记
            // 这样下一次 CursorPosCallback 就不会计算错误的 Delta
            ResetMouse();
            
            // 重新锁定光标 (有些系统切出去会自动解锁光标)
            glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
        } else {
            // 窗口失去了焦点
            // 可以选择解锁光标，让玩家能操作其他窗口
            glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
        }
    }

    void Input::ResetMouse() {
        s_firstMouseInput = true;
        
        // 可选：顺便清空一下 Delta，防止残留
        m_mouseDelta = Vector2(0, 0); 
    }

    //获取数据信息

    //传递键盘信息
    bool Input::GetKey(Key key){
        return m_keys[(int)key];
    }

    bool Input::GetKeyDown(Key key){
        return m_keys[(int)key] && !m_lastKeys[(int)key];
    }

    bool Input::GetKeyUp(Key key){
        return !m_keys[(int)key] && m_lastKeys[(int)key];
    }

    //传递鼠标信息
    bool Input::GetMouseButton(MouseButton button){
        return m_mouseButtons[(int)button];
    }

    bool Input::GetMouseButtonDown(MouseButton button){
        return m_mouseButtons[(int)button] && !m_lastMouseButtons[(int)button];
    }

    bool Input::GetMouseButtonUp(MouseButton button){
        return !m_mouseButtons[(int)button] && m_lastMouseButtons[(int)button];
    }

    //传递参数
    Vector2 Input::GetMousePosition(){
        return m_mousePos;
    }

    float Input::GetMouseScroll(){
        return m_scrollDelta;
    }

    Vector2 Input::GetMouseDelta(){
        return m_mouseDelta;
    }
}
/*
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
*/
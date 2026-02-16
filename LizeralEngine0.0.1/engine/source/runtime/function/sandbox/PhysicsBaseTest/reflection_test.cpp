#include <iostream>
#include <string>
#include <cassert>
#include <vector>

// 1. GLAD/GLFW (必须放在最前面)
#include <glad/glad.h>
#include <GLFW/glfw3.h>

// 2. 引入引擎基础组件
#include "runtime/function/ecs/components/Transform/TransformComponent.h"
#include "runtime/core/math/vector3.h"

// =======================================================
// 【Mock 区域】模拟序列化器
// 目的：为了单独测试反射逻辑，避免编译庞大的 Serializer 实现代码
// =======================================================

// 1. 引入真正的 PJson 定义 (避免重定义错误)
#include "runtime/core/meta/json.h" 

namespace Lizeral {
    // 2. 定义 PSerializer
    // 虽然 component.h 里可能前置声明了它，但我们没有链接 serializer.cpp
    // 所以在这里提供一个空的内联实现，骗过链接器
    class PSerializer {
    public:
        template<typename T> 
        static void read(const PJson&, T&) {}
        
        template<typename T> 
        static PJson write(const T&) { return PJson(); }
    };
}
// =======================================================

// 3. 引入生成的反射代码
// (此时 PSerializer 已定义，PJson 已引入，编译器不会报错)
#include "_generated/reflection/TransformComponent.reflection.gen.h" 

// 4. 其他引擎头文件
#include "runtime/function/ecs/registry.h"
#include "runtime/function/physics/PhysicsSystem.h"
#include "runtime/function/physics/PhysicsScene.h"

using namespace Lizeral;

// ========================================================================
// 模拟编辑器行为的测试函数
// ========================================================================
void TestTransformReflection() {
    std::cout << "========== 开始反射系统测试 ==========" << std::endl;

    // 1. 初始化反射注册
    // 在引擎正式启动时，会自动调用所有生成的注册函数
    // 这里手动调用一次以确保环境就绪
    std::cout << "[Step 1] 注册反射元数据..." << std::endl;
    Lizeral::Reflection::TypeWrappersRegister::TransformComponent();

    // 2. 创建一个运行时实例
    std::cout << "[Step 2] 创建 TransformComponent 实例..." << std::endl;
    auto* transformComp = new TransformComponent();
    
    // 初始化默认值
    transformComp->setPosition(Vector3(0, 0, 0));
    // 注意：如果有 cleanDirty() 接口，建议在此处调用以重置脏标记

    // 3. [核心测试] 模拟编辑器：通过字符串 "m_position" 修改值
    std::cout << "[Step 3] 模拟编辑器修改属性..." << std::endl;
    
    // 假设这是编辑器面板输入的新数值
    Vector3 newEditorPosition(100.0f, 50.0f, -20.0f);
    
    // A. 获取反射 Operator (模拟全局查找)
    using Operator = Lizeral::Reflection::TypeFieldReflectionOparator::TypeTransformComponentOperator;

    // B. 验证字段名是否存在
    std::string targetFieldName = "m_position";
    if (std::string(Operator::getFieldName_m_position()) == targetFieldName) {
        std::cout << "  -> 成功找到字段: " << targetFieldName << std::endl;
        std::cout << "  -> 字段类型: " << Operator::getFieldTypeName_m_position() << std::endl;

        // C. 通过反射进行 Set 操作
        // 编辑器通常持有 void* 指针
        void* instancePtr = static_cast<void*>(transformComp);
        void* valuePtr = static_cast<void*>(&newEditorPosition);

        // 调用生成的 set 函数 -> 内部调用 TryNotifyReflectionUpdated -> 触发 setDirty
        Operator::set_m_position(instancePtr, valuePtr);
        
        std::cout << "  -> 已通过反射 Set 设置新值。" << std::endl;
    }

    // 4. [验证结果] 检查数值是否改变
    std::cout << "[Step 4] 验证数值一致性..." << std::endl;
    Vector3 currentPos = transformComp->getPosition();
    bool valMatch = (currentPos.x == 100.0f && currentPos.y == 50.0f && currentPos.z == -20.0f);
    
    if (valMatch) {
        std::cout << "  [PASS] 数值更新成功！Transform 现在的坐标是: " 
                  << currentPos.x << ", " << currentPos.y << ", " << currentPos.z << std::endl;
    } else {
        std::cout << "  [FAIL] 数值更新失败！" << std::endl;
    }

    // 5. [关键验证] 检查脏标记 (Dirty Flag) 是否被触发
    std::cout << "[Step 5] 验证脏标记联动 (Reflection -> DirtyFlag)..." << std::endl;
    
    bool isDirty = transformComp->isDirty(Lizeral::TRANS_DIRTY_POSITION);
    
    if (isDirty) {
        std::cout << "  [PASS] 脏标记已触发！(TRANS_DIRTY_POSITION)" << std::endl;
        std::cout << "  -> 物理系统将在下一帧同步此数据。" << std::endl;
    } else {
        std::cout << "  [FAIL] 脏标记未触发！请检查 REFLECTION_BIND_DIRTY 宏定义。" << std::endl;
    }

    // 6. 清理
    delete transformComp;
    std::cout << "========== 测试结束 ==========" << std::endl;
}

int main() {
    TestTransformReflection();
    return 0;
}
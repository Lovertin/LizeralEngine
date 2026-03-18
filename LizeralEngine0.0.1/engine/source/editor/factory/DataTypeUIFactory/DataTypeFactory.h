#pragma once
#include <QWidget>
#include <unordered_map>
#include <string>
#include <memory>
// 请确保这里包含了你的反射系统头文件，因为我们需要用到 FieldAccessor
#include "runtime/core/meta/reflection/reflection.h"
#include "editor/selection/EditorSelection.h"

#include "editor/context/EditorContext.h"
#include "editor/command/EditPropertyCommand.h"

namespace Lizeral {

    // --- 1. 基础绘制器接口 ---
    class IPropertyDrawer {
    public:
        virtual ~IPropertyDrawer() = default;

        // 核心接口：负责生成 UI 并绑定数据
        virtual QWidget* DrawProperty(Reflection::FieldAccessor& accessor, void* instance, QWidget* parent) = 0;
    };

    // --- 2. 数据类型 UI 工厂 ---
    class DataTypeUIFactory {
    public:
        // 在编辑器启动时（或第一次需要时）调用，注册所有微型工厂
        static void Initialize();

        // 核心路由函数：根据 accessor 的类型和标签，分发给具体的 Drawer
        static QWidget* CreateFieldWidget(Reflection::FieldAccessor& accessor, void* instance, QWidget* parent = nullptr);

    private:
    
        static std::unordered_map<std::string, std::unique_ptr<IPropertyDrawer>> s_DrawerRegistry;
        static bool s_Initialized;
    };

} // namespace Lizeral
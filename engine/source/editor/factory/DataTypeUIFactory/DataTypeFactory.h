#pragma once
#include <QWidget>
#include <unordered_map>
#include <string>
#include <memory>
#include "runtime/core/meta/reflection/reflection.h"
#include "editor/selection/EditorSelection.h"

#include "editor/context/EditorContext.h"
#include "editor/command/EditPropertyCommand.h"

namespace Lizeral {

    // basic drawer interface
    class IPropertyDrawer {
    public:
        virtual ~IPropertyDrawer() = default;

        virtual QWidget* DrawProperty(Reflection::FieldAccessor& accessor, void* instance, QWidget* parent) = 0;
    };

    class DataTypeUIFactory {
    public:
        static void Initialize();

        static QWidget* CreateFieldWidget(Reflection::FieldAccessor& accessor, void* instance, QWidget* parent = nullptr);

    private:
    
        static std::unordered_map<std::string, std::unique_ptr<IPropertyDrawer>> s_DrawerRegistry;
        static bool s_Initialized;
    };

} // namespace Lizeral
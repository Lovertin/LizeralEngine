// ComponentUIFactory.cpp
#include "ComponentUIFactory.h"

namespace Lizeral{

    QWidget* ComponentUIFactory::CreateWidgetFor(Component* component, QWidget* parent) {
        if (!component) return nullptr;

        // 伪代码：这里通常通过反射的类型名，或者 RTTI (dynamic_cast) 来判断组件类型
        // std::string typeName = component->GetTypeName();
        
        // 假设我们识别出这是一个 TransformComponent
        // if (typeName == "TransformComponent") {
        //     return CreateTransformUI(static_cast<TransformComponent*>(component), parent);
        // }

        // 兜底 UI：如果我们不认识这个组件，只显示它的名字
        QGroupBox* box = new QGroupBox("Unknown Component", parent);
        QVBoxLayout* layout = new QVBoxLayout(box);
        layout->addWidget(new QLabel("No UI builder registered for this component.", box));
        return box;
    }
}
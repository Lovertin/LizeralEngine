// ComponentUIFactory.h
#pragma once
#include <QWidget>
#include <QVBoxLayout>
#include <QLabel>
#include <QGroupBox>

#include "runtime/function/ecs/components/component.h"
#include "runtime/function/ecs/components/componentAll.h"

namespace Lizeral{

    class ComponentUIFactory {
    public:

        static QWidget* CreateWidgetFor(Component* component, QWidget* parent = nullptr);

    private:

        static QWidget* CreateTransformUI(TransformComponent* transform, QWidget* parent);
    };

}
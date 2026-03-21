// ComponentUIFactory.h
#pragma once
#include <QWidget>
#include "runtime/core/meta/reflection/reflection.h"


namespace Lizeral{

    class ComponentUIFactory {
    public:
        static QWidget* CreateComponentWidget(Reflection::TypeMeta& meta,void* instance,QWidget* parent = nullptr);
    };

}
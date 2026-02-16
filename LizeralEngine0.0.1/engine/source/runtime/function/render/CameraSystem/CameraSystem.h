#pragma once
#include "runtime/function/ecs/registry.h"

namespace Lizeral{
    class CameraSystem{

    public:
        CameraSystem();
        virtual ~CameraSystem();

        void Tick(Registry& registry);

    private:
        

    };
}
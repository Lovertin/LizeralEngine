#pragma once
#include "runtime/core/meta/reflection/reflection.h"

namespace Lizeral{
    REFLECTION_TYPE(Test);
    CLASS(Test,Fields){
        REFLECTION_BODY(Test)

    public:
        float x,y;
        int z;
        std::string s;
        float getX(){ return x;}
        float getY(){ return y; }
    };
}
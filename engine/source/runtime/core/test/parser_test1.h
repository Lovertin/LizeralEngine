#pragma once
#include "runtime/core/meta/reflection/reflection.h"

class Test1{
public:
    float x,y;
    Test1():x(0.0f),y(0.0f){}

    void setX(int val){
        x=val;
    }
};
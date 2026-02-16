#pragma once
#include "runtime/core/meta/reflection/reflection.h"
#include "_generated/serializer/all_serializer.h"
#include "_generated\reflection\TransformComponent.reflection.gen.h"
#include "_generated\reflection\quaternion.reflection.gen.h"
#include "_generated\reflection\vector3.reflection.gen.h"
#include "_generated\reflection\color.reflection.gen.h"
#include "_generated\reflection\vector4.reflection.gen.h"
#include "_generated\reflection\axis_aligned.reflection.gen.h"
#include "_generated\reflection\matrix4.reflection.gen.h"
#include "_generated\reflection\vector2.reflection.gen.h"
#include "_generated\reflection\parser_test.reflection.gen.h"
#include "_generated\reflection\component.reflection.gen.h"
#include "_generated\reflection\transform.reflection.gen.h"
#include "_generated\reflection\ColliderComponent.reflection.gen.h"
#include "_generated\reflection\CameraComponent.reflection.gen.h"
#include "_generated\reflection\RigidBodyComponent.reflection.gen.h"

namespace Lizeral{
namespace Reflection{
    void TypeMetaRegister::Register(){
    TypeWrappersRegister::TransformComponent();
    TypeWrappersRegister::Quaternion();
    TypeWrappersRegister::Vector3();
    TypeWrappersRegister::Color();
    TypeWrappersRegister::Vector4();
    TypeWrappersRegister::AxisAlignedBox();
    TypeWrappersRegister::Matrix4x4_();
    TypeWrappersRegister::Vector2();
    TypeWrappersRegister::Test();
    TypeWrappersRegister::Component();
    TypeWrappersRegister::Transform();
    TypeWrappersRegister::ColliderComponent();
    TypeWrappersRegister::CameraComponent();
    TypeWrappersRegister::RigidBodyComponent();
    }
}
}


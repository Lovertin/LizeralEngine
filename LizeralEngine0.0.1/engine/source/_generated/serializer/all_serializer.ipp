#pragma once
#include "_generated\serializer\PBRMaterial.serializer.gen.h"
#include "_generated\serializer\Texture2D.serializer.gen.h"
#include "_generated\serializer\TransformComponent.serializer.gen.h"
#include "_generated\serializer\RigidBodyComponent.serializer.gen.h"
#include "_generated\serializer\TextureCube.serializer.gen.h"
#include "_generated\serializer\ModelComponent.serializer.gen.h"
#include "_generated\serializer\resource.serializer.gen.h"
#include "_generated\serializer\Mesh.serializer.gen.h"
#include "_generated\serializer\quaternion.serializer.gen.h"
#include "_generated\serializer\DirectionalLightComponent.serializer.gen.h"
#include "_generated\serializer\vector3.serializer.gen.h"
#include "_generated\serializer\color.serializer.gen.h"
#include "_generated\serializer\vector4.serializer.gen.h"
#include "_generated\serializer\axis_aligned.serializer.gen.h"
#include "_generated\serializer\CameraControlComponent.serializer.gen.h"
#include "_generated\serializer\matrix4.serializer.gen.h"
#include "_generated\serializer\vector2.serializer.gen.h"
#include "_generated\serializer\parser_test.serializer.gen.h"
#include "_generated\serializer\component.serializer.gen.h"
#include "_generated\serializer\transform.serializer.gen.h"
#include "_generated\serializer\CameraComponent.serializer.gen.h"
#include "_generated\serializer\ColliderComponent.serializer.gen.h"
namespace Lizeral{
    template<>
    PJson PSerializer::write(const PBRMaterial& instance){
        PJson::object  ret_context;
        auto&&  json_context_0 = PSerializer::write(*(Lizeral::Material*)&instance);
        assert(json_context_0.is_object());
        auto&& json_context_map_0 = json_context_0.object_items();
        ret_context.insert(json_context_map_0.begin() , json_context_map_0.end());
        
        return  PJson(ret_context);
    }
    template<>
    PBRMaterial& PSerializer::read(const PJson& json_context, PBRMaterial& instance){
        assert(json_context.is_object());
        PSerializer::read(json_context,*(Lizeral::Material*)&instance);
        
        return instance;
    }
    template<>
    PJson PSerializer::write(const Texture2D& instance){
        PJson::object  ret_context;
        auto&&  json_context_0 = PSerializer::write(*(Lizeral::Resource*)&instance);
        assert(json_context_0.is_object());
        auto&& json_context_map_0 = json_context_0.object_items();
        ret_context.insert(json_context_map_0.begin() , json_context_map_0.end());
        
        return  PJson(ret_context);
    }
    template<>
    Texture2D& PSerializer::read(const PJson& json_context, Texture2D& instance){
        assert(json_context.is_object());
        PSerializer::read(json_context,*(Lizeral::Resource*)&instance);
        
        return instance;
    }
    template<>
    PJson PSerializer::write(const TransformComponent& instance){
        PJson::object  ret_context;
        auto&&  json_context_0 = PSerializer::write(*(Lizeral::Component*)&instance);
        assert(json_context_0.is_object());
        auto&& json_context_map_0 = json_context_0.object_items();
        ret_context.insert(json_context_map_0.begin() , json_context_map_0.end());
        ret_context.insert_or_assign("position", PSerializer::write(instance.m_position));
        ret_context.insert_or_assign("rotation", PSerializer::write(instance.m_rotation));
        ret_context.insert_or_assign("scale", PSerializer::write(instance.m_scale));
        return  PJson(ret_context);
    }
    template<>
    TransformComponent& PSerializer::read(const PJson& json_context, TransformComponent& instance){
        assert(json_context.is_object());
        PSerializer::read(json_context,*(Lizeral::Component*)&instance);
        if(!json_context["position"].is_null()){
            PSerializer::read(json_context["position"], instance.m_position);
        }
        if(!json_context["rotation"].is_null()){
            PSerializer::read(json_context["rotation"], instance.m_rotation);
        }
        if(!json_context["scale"].is_null()){
            PSerializer::read(json_context["scale"], instance.m_scale);
        }
        return instance;
    }
    template<>
    PJson PSerializer::write(const RigidBodyComponent& instance){
        PJson::object  ret_context;
        auto&&  json_context_0 = PSerializer::write(*(Lizeral::Component*)&instance);
        assert(json_context_0.is_object());
        auto&& json_context_map_0 = json_context_0.object_items();
        ret_context.insert(json_context_map_0.begin() , json_context_map_0.end());
        ret_context.insert_or_assign("mass", PSerializer::write(instance.m_mass));
        ret_context.insert_or_assign("friction", PSerializer::write(instance.m_friction));
        ret_context.insert_or_assign("is_kinematic", PSerializer::write(instance.m_is_kinematic));
        ret_context.insert_or_assign("restitution", PSerializer::write(instance.m_restitution));
        return  PJson(ret_context);
    }
    template<>
    RigidBodyComponent& PSerializer::read(const PJson& json_context, RigidBodyComponent& instance){
        assert(json_context.is_object());
        PSerializer::read(json_context,*(Lizeral::Component*)&instance);
        if(!json_context["mass"].is_null()){
            PSerializer::read(json_context["mass"], instance.m_mass);
        }
        if(!json_context["friction"].is_null()){
            PSerializer::read(json_context["friction"], instance.m_friction);
        }
        if(!json_context["is_kinematic"].is_null()){
            PSerializer::read(json_context["is_kinematic"], instance.m_is_kinematic);
        }
        if(!json_context["restitution"].is_null()){
            PSerializer::read(json_context["restitution"], instance.m_restitution);
        }
        return instance;
    }
    template<>
    PJson PSerializer::write(const TextureCube& instance){
        PJson::object  ret_context;
        auto&&  json_context_0 = PSerializer::write(*(Lizeral::Resource*)&instance);
        assert(json_context_0.is_object());
        auto&& json_context_map_0 = json_context_0.object_items();
        ret_context.insert(json_context_map_0.begin() , json_context_map_0.end());
        
        return  PJson(ret_context);
    }
    template<>
    TextureCube& PSerializer::read(const PJson& json_context, TextureCube& instance){
        assert(json_context.is_object());
        PSerializer::read(json_context,*(Lizeral::Resource*)&instance);
        
        return instance;
    }
    template<>
    PJson PSerializer::write(const ModelComponent& instance){
        PJson::object  ret_context;
        auto&&  json_context_0 = PSerializer::write(*(Lizeral::Component*)&instance);
        assert(json_context_0.is_object());
        auto&& json_context_map_0 = json_context_0.object_items();
        ret_context.insert(json_context_map_0.begin() , json_context_map_0.end());
        
        return  PJson(ret_context);
    }
    template<>
    ModelComponent& PSerializer::read(const PJson& json_context, ModelComponent& instance){
        assert(json_context.is_object());
        PSerializer::read(json_context,*(Lizeral::Component*)&instance);
        
        return instance;
    }
    template<>
    PJson PSerializer::write(const Resource& instance){
        PJson::object  ret_context;
        
        ret_context.insert_or_assign("path", PSerializer::write(instance.m_path));
        return  PJson(ret_context);
    }
    template<>
    Resource& PSerializer::read(const PJson& json_context, Resource& instance){
        assert(json_context.is_object());
        
        if(!json_context["path"].is_null()){
            PSerializer::read(json_context["path"], instance.m_path);
        }
        return instance;
    }
    template<>
    PJson PSerializer::write(const Mesh& instance){
        PJson::object  ret_context;
        
        
        return  PJson(ret_context);
    }
    template<>
    Mesh& PSerializer::read(const PJson& json_context, Mesh& instance){
        assert(json_context.is_object());
        
        
        return instance;
    }
    template<>
    PJson PSerializer::write(const Quaternion& instance){
        PJson::object  ret_context;
        
        ret_context.insert_or_assign("w", PSerializer::write(instance.w));
        ret_context.insert_or_assign("x", PSerializer::write(instance.x));
        ret_context.insert_or_assign("y", PSerializer::write(instance.y));
        ret_context.insert_or_assign("z", PSerializer::write(instance.z));
        return  PJson(ret_context);
    }
    template<>
    Quaternion& PSerializer::read(const PJson& json_context, Quaternion& instance){
        assert(json_context.is_object());
        
        if(!json_context["w"].is_null()){
            PSerializer::read(json_context["w"], instance.w);
        }
        if(!json_context["x"].is_null()){
            PSerializer::read(json_context["x"], instance.x);
        }
        if(!json_context["y"].is_null()){
            PSerializer::read(json_context["y"], instance.y);
        }
        if(!json_context["z"].is_null()){
            PSerializer::read(json_context["z"], instance.z);
        }
        return instance;
    }
    template<>
    PJson PSerializer::write(const DirectionLightComponent& instance){
        PJson::object  ret_context;
        auto&&  json_context_0 = PSerializer::write(*(Lizeral::Component*)&instance);
        assert(json_context_0.is_object());
        auto&& json_context_map_0 = json_context_0.object_items();
        ret_context.insert(json_context_map_0.begin() , json_context_map_0.end());
        ret_context.insert_or_assign("color", PSerializer::write(instance.m_color));
        ret_context.insert_or_assign("intensity", PSerializer::write(instance.m_intensity));
        ret_context.insert_or_assign("isGlobal", PSerializer::write(instance.m_isGlobal));
        return  PJson(ret_context);
    }
    template<>
    DirectionLightComponent& PSerializer::read(const PJson& json_context, DirectionLightComponent& instance){
        assert(json_context.is_object());
        PSerializer::read(json_context,*(Lizeral::Component*)&instance);
        if(!json_context["color"].is_null()){
            PSerializer::read(json_context["color"], instance.m_color);
        }
        if(!json_context["intensity"].is_null()){
            PSerializer::read(json_context["intensity"], instance.m_intensity);
        }
        if(!json_context["isGlobal"].is_null()){
            PSerializer::read(json_context["isGlobal"], instance.m_isGlobal);
        }
        return instance;
    }
    template<>
    PJson PSerializer::write(const Vector3& instance){
        PJson::object  ret_context;
        
        ret_context.insert_or_assign("x", PSerializer::write(instance.x));
        ret_context.insert_or_assign("y", PSerializer::write(instance.y));
        ret_context.insert_or_assign("z", PSerializer::write(instance.z));
        return  PJson(ret_context);
    }
    template<>
    Vector3& PSerializer::read(const PJson& json_context, Vector3& instance){
        assert(json_context.is_object());
        
        if(!json_context["x"].is_null()){
            PSerializer::read(json_context["x"], instance.x);
        }
        if(!json_context["y"].is_null()){
            PSerializer::read(json_context["y"], instance.y);
        }
        if(!json_context["z"].is_null()){
            PSerializer::read(json_context["z"], instance.z);
        }
        return instance;
    }
    template<>
    PJson PSerializer::write(const Color& instance){
        PJson::object  ret_context;
        
        ret_context.insert_or_assign("r", PSerializer::write(instance.r));
        ret_context.insert_or_assign("g", PSerializer::write(instance.g));
        ret_context.insert_or_assign("b", PSerializer::write(instance.b));
        return  PJson(ret_context);
    }
    template<>
    Color& PSerializer::read(const PJson& json_context, Color& instance){
        assert(json_context.is_object());
        
        if(!json_context["r"].is_null()){
            PSerializer::read(json_context["r"], instance.r);
        }
        if(!json_context["g"].is_null()){
            PSerializer::read(json_context["g"], instance.g);
        }
        if(!json_context["b"].is_null()){
            PSerializer::read(json_context["b"], instance.b);
        }
        return instance;
    }
    template<>
    PJson PSerializer::write(const Vector4& instance){
        PJson::object  ret_context;
        
        ret_context.insert_or_assign("x", PSerializer::write(instance.x));
        ret_context.insert_or_assign("y", PSerializer::write(instance.y));
        ret_context.insert_or_assign("z", PSerializer::write(instance.z));
        ret_context.insert_or_assign("w", PSerializer::write(instance.w));
        return  PJson(ret_context);
    }
    template<>
    Vector4& PSerializer::read(const PJson& json_context, Vector4& instance){
        assert(json_context.is_object());
        
        if(!json_context["x"].is_null()){
            PSerializer::read(json_context["x"], instance.x);
        }
        if(!json_context["y"].is_null()){
            PSerializer::read(json_context["y"], instance.y);
        }
        if(!json_context["z"].is_null()){
            PSerializer::read(json_context["z"], instance.z);
        }
        if(!json_context["w"].is_null()){
            PSerializer::read(json_context["w"], instance.w);
        }
        return instance;
    }
    template<>
    PJson PSerializer::write(const AxisAlignedBox& instance){
        PJson::object  ret_context;
        
        ret_context.insert_or_assign("center", PSerializer::write(instance.m_center));
        ret_context.insert_or_assign("half_extent", PSerializer::write(instance.m_half_extent));
        ret_context.insert_or_assign("min_corner", PSerializer::write(instance.m_min_corner));
        ret_context.insert_or_assign("max_corner", PSerializer::write(instance.m_max_corner));
        return  PJson(ret_context);
    }
    template<>
    AxisAlignedBox& PSerializer::read(const PJson& json_context, AxisAlignedBox& instance){
        assert(json_context.is_object());
        
        if(!json_context["center"].is_null()){
            PSerializer::read(json_context["center"], instance.m_center);
        }
        if(!json_context["half_extent"].is_null()){
            PSerializer::read(json_context["half_extent"], instance.m_half_extent);
        }
        if(!json_context["min_corner"].is_null()){
            PSerializer::read(json_context["min_corner"], instance.m_min_corner);
        }
        if(!json_context["max_corner"].is_null()){
            PSerializer::read(json_context["max_corner"], instance.m_max_corner);
        }
        return instance;
    }
    template<>
    PJson PSerializer::write(const CameraControlComponent& instance){
        PJson::object  ret_context;
        auto&&  json_context_0 = PSerializer::write(*(Lizeral::Component*)&instance);
        assert(json_context_0.is_object());
        auto&& json_context_map_0 = json_context_0.object_items();
        ret_context.insert(json_context_map_0.begin() , json_context_map_0.end());
        ret_context.insert_or_assign("move_speed", PSerializer::write(instance.move_speed));
        ret_context.insert_or_assign("speedMultiplier", PSerializer::write(instance.m_speedMultiplier));
        ret_context.insert_or_assign("sensitivityX", PSerializer::write(instance.m_sensitivityX));
        ret_context.insert_or_assign("sensitivityY", PSerializer::write(instance.m_sensitivityY));
        ret_context.insert_or_assign("pitch", PSerializer::write(instance.m_pitch));
        ret_context.insert_or_assign("yaw", PSerializer::write(instance.m_yaw));
        return  PJson(ret_context);
    }
    template<>
    CameraControlComponent& PSerializer::read(const PJson& json_context, CameraControlComponent& instance){
        assert(json_context.is_object());
        PSerializer::read(json_context,*(Lizeral::Component*)&instance);
        if(!json_context["move_speed"].is_null()){
            PSerializer::read(json_context["move_speed"], instance.move_speed);
        }
        if(!json_context["speedMultiplier"].is_null()){
            PSerializer::read(json_context["speedMultiplier"], instance.m_speedMultiplier);
        }
        if(!json_context["sensitivityX"].is_null()){
            PSerializer::read(json_context["sensitivityX"], instance.m_sensitivityX);
        }
        if(!json_context["sensitivityY"].is_null()){
            PSerializer::read(json_context["sensitivityY"], instance.m_sensitivityY);
        }
        if(!json_context["pitch"].is_null()){
            PSerializer::read(json_context["pitch"], instance.m_pitch);
        }
        if(!json_context["yaw"].is_null()){
            PSerializer::read(json_context["yaw"], instance.m_yaw);
        }
        return instance;
    }
    template<>
    PJson PSerializer::write(const Matrix4x4_& instance){
        PJson::object  ret_context;
        
        ret_context.insert_or_assign("v0", PSerializer::write(instance.v0));
        ret_context.insert_or_assign("v1", PSerializer::write(instance.v1));
        ret_context.insert_or_assign("v2", PSerializer::write(instance.v2));
        ret_context.insert_or_assign("v3", PSerializer::write(instance.v3));
        ret_context.insert_or_assign("v4", PSerializer::write(instance.v4));
        ret_context.insert_or_assign("v5", PSerializer::write(instance.v5));
        ret_context.insert_or_assign("v6", PSerializer::write(instance.v6));
        ret_context.insert_or_assign("v7", PSerializer::write(instance.v7));
        ret_context.insert_or_assign("v8", PSerializer::write(instance.v8));
        ret_context.insert_or_assign("v9", PSerializer::write(instance.v9));
        ret_context.insert_or_assign("v10", PSerializer::write(instance.v10));
        ret_context.insert_or_assign("v11", PSerializer::write(instance.v11));
        ret_context.insert_or_assign("v12", PSerializer::write(instance.v12));
        ret_context.insert_or_assign("v13", PSerializer::write(instance.v13));
        ret_context.insert_or_assign("v14", PSerializer::write(instance.v14));
        ret_context.insert_or_assign("v15", PSerializer::write(instance.v15));
        return  PJson(ret_context);
    }
    template<>
    Matrix4x4_& PSerializer::read(const PJson& json_context, Matrix4x4_& instance){
        assert(json_context.is_object());
        
        if(!json_context["v0"].is_null()){
            PSerializer::read(json_context["v0"], instance.v0);
        }
        if(!json_context["v1"].is_null()){
            PSerializer::read(json_context["v1"], instance.v1);
        }
        if(!json_context["v2"].is_null()){
            PSerializer::read(json_context["v2"], instance.v2);
        }
        if(!json_context["v3"].is_null()){
            PSerializer::read(json_context["v3"], instance.v3);
        }
        if(!json_context["v4"].is_null()){
            PSerializer::read(json_context["v4"], instance.v4);
        }
        if(!json_context["v5"].is_null()){
            PSerializer::read(json_context["v5"], instance.v5);
        }
        if(!json_context["v6"].is_null()){
            PSerializer::read(json_context["v6"], instance.v6);
        }
        if(!json_context["v7"].is_null()){
            PSerializer::read(json_context["v7"], instance.v7);
        }
        if(!json_context["v8"].is_null()){
            PSerializer::read(json_context["v8"], instance.v8);
        }
        if(!json_context["v9"].is_null()){
            PSerializer::read(json_context["v9"], instance.v9);
        }
        if(!json_context["v10"].is_null()){
            PSerializer::read(json_context["v10"], instance.v10);
        }
        if(!json_context["v11"].is_null()){
            PSerializer::read(json_context["v11"], instance.v11);
        }
        if(!json_context["v12"].is_null()){
            PSerializer::read(json_context["v12"], instance.v12);
        }
        if(!json_context["v13"].is_null()){
            PSerializer::read(json_context["v13"], instance.v13);
        }
        if(!json_context["v14"].is_null()){
            PSerializer::read(json_context["v14"], instance.v14);
        }
        if(!json_context["v15"].is_null()){
            PSerializer::read(json_context["v15"], instance.v15);
        }
        return instance;
    }
    template<>
    PJson PSerializer::write(const Vector2& instance){
        PJson::object  ret_context;
        
        ret_context.insert_or_assign("x", PSerializer::write(instance.x));
        ret_context.insert_or_assign("y", PSerializer::write(instance.y));
        return  PJson(ret_context);
    }
    template<>
    Vector2& PSerializer::read(const PJson& json_context, Vector2& instance){
        assert(json_context.is_object());
        
        if(!json_context["x"].is_null()){
            PSerializer::read(json_context["x"], instance.x);
        }
        if(!json_context["y"].is_null()){
            PSerializer::read(json_context["y"], instance.y);
        }
        return instance;
    }
    template<>
    PJson PSerializer::write(const Test& instance){
        PJson::object  ret_context;
        
        ret_context.insert_or_assign("x", PSerializer::write(instance.x));
        ret_context.insert_or_assign("y", PSerializer::write(instance.y));
        ret_context.insert_or_assign("z", PSerializer::write(instance.z));
        ret_context.insert_or_assign("s", PSerializer::write(instance.s));
        return  PJson(ret_context);
    }
    template<>
    Test& PSerializer::read(const PJson& json_context, Test& instance){
        assert(json_context.is_object());
        
        if(!json_context["x"].is_null()){
            PSerializer::read(json_context["x"], instance.x);
        }
        if(!json_context["y"].is_null()){
            PSerializer::read(json_context["y"], instance.y);
        }
        if(!json_context["z"].is_null()){
            PSerializer::read(json_context["z"], instance.z);
        }
        if(!json_context["s"].is_null()){
            PSerializer::read(json_context["s"], instance.s);
        }
        return instance;
    }
    template<>
    PJson PSerializer::write(const Component& instance){
        PJson::object  ret_context;
        
        
        return  PJson(ret_context);
    }
    template<>
    Component& PSerializer::read(const PJson& json_context, Component& instance){
        assert(json_context.is_object());
        
        
        return instance;
    }
    template<>
    PJson PSerializer::write(const Transform& instance){
        PJson::object  ret_context;
        
        ret_context.insert_or_assign("position", PSerializer::write(instance.m_position));
        ret_context.insert_or_assign("scale", PSerializer::write(instance.m_scale));
        ret_context.insert_or_assign("rotation", PSerializer::write(instance.m_rotation));
        return  PJson(ret_context);
    }
    template<>
    Transform& PSerializer::read(const PJson& json_context, Transform& instance){
        assert(json_context.is_object());
        
        if(!json_context["position"].is_null()){
            PSerializer::read(json_context["position"], instance.m_position);
        }
        if(!json_context["scale"].is_null()){
            PSerializer::read(json_context["scale"], instance.m_scale);
        }
        if(!json_context["rotation"].is_null()){
            PSerializer::read(json_context["rotation"], instance.m_rotation);
        }
        return instance;
    }
    template<>
    PJson PSerializer::write(const CameraComponent& instance){
        PJson::object  ret_context;
        auto&&  json_context_0 = PSerializer::write(*(Lizeral::Component*)&instance);
        assert(json_context_0.is_object());
        auto&& json_context_map_0 = json_context_0.object_items();
        ret_context.insert(json_context_map_0.begin() , json_context_map_0.end());
        ret_context.insert_or_assign("fov", PSerializer::write(instance.m_fov));
        ret_context.insert_or_assign("aspect", PSerializer::write(instance.m_aspect));
        ret_context.insert_or_assign("zNear", PSerializer::write(instance.m_zNear));
        ret_context.insert_or_assign("zFar", PSerializer::write(instance.m_zFar));
        ret_context.insert_or_assign("isMainCamera", PSerializer::write(instance.isMainCamera));
        return  PJson(ret_context);
    }
    template<>
    CameraComponent& PSerializer::read(const PJson& json_context, CameraComponent& instance){
        assert(json_context.is_object());
        PSerializer::read(json_context,*(Lizeral::Component*)&instance);
        if(!json_context["fov"].is_null()){
            PSerializer::read(json_context["fov"], instance.m_fov);
        }
        if(!json_context["aspect"].is_null()){
            PSerializer::read(json_context["aspect"], instance.m_aspect);
        }
        if(!json_context["zNear"].is_null()){
            PSerializer::read(json_context["zNear"], instance.m_zNear);
        }
        if(!json_context["zFar"].is_null()){
            PSerializer::read(json_context["zFar"], instance.m_zFar);
        }
        if(!json_context["isMainCamera"].is_null()){
            PSerializer::read(json_context["isMainCamera"], instance.isMainCamera);
        }
        return instance;
    }
    template<>
    PJson PSerializer::write(const ColliderComponent& instance){
        PJson::object  ret_context;
        auto&&  json_context_0 = PSerializer::write(*(Lizeral::Component*)&instance);
        assert(json_context_0.is_object());
        auto&& json_context_map_0 = json_context_0.object_items();
        ret_context.insert(json_context_map_0.begin() , json_context_map_0.end());
        ret_context.insert_or_assign("type", PSerializer::write(instance.m_type));
        ret_context.insert_or_assign("size", PSerializer::write(instance.m_size));
        ret_context.insert_or_assign("radius", PSerializer::write(instance.m_radius));
        ret_context.insert_or_assign("height", PSerializer::write(instance.m_height));
        ret_context.insert_or_assign("offset", PSerializer::write(instance.m_offset));
        return  PJson(ret_context);
    }
    template<>
    ColliderComponent& PSerializer::read(const PJson& json_context, ColliderComponent& instance){
        assert(json_context.is_object());
        PSerializer::read(json_context,*(Lizeral::Component*)&instance);
        if(!json_context["type"].is_null()){
            PSerializer::read(json_context["type"], instance.m_type);
        }
        if(!json_context["size"].is_null()){
            PSerializer::read(json_context["size"], instance.m_size);
        }
        if(!json_context["radius"].is_null()){
            PSerializer::read(json_context["radius"], instance.m_radius);
        }
        if(!json_context["height"].is_null()){
            PSerializer::read(json_context["height"], instance.m_height);
        }
        if(!json_context["offset"].is_null()){
            PSerializer::read(json_context["offset"], instance.m_offset);
        }
        return instance;
    }

}


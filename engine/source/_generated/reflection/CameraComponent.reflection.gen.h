#pragma once
#include "runtime/core/meta/serializer/serializer.h"
#include "runtime\function\ecs\components\Camera\CameraComponent.h"

namespace Lizeral{
    class CameraComponent;
namespace Reflection{
namespace TypeFieldReflectionOparator{
    class TypeCameraComponentOperator{
    public:
        static const char* getClassName(){ return "CameraComponent";}
        static void* constructorWithJson(const PJson& json_context){
            CameraComponent* ret_instance= new CameraComponent;
            PSerializer::read(json_context, *ret_instance);
            return ret_instance;
        }
        static PJson writeByName(void* instance){
            return PSerializer::write(*(CameraComponent*)instance);
        }
        // base class
        static std::vector<ReflectionInstance> getCameraComponentBaseClassReflectionInstanceList(void* instance){
            std::vector<ReflectionInstance> out_list;
            out_list.push_back(TypeMetaDef(Lizeral::Component, static_cast<CameraComponent*>(instance)));
            return out_list;
        }
        // fields
        static const char* getFieldName_m_fov(){ return "m_fov";}
        static const char* getFieldTypeName_m_fov(){ return "float";}
        static void set_m_fov(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<CameraComponent*>(instance);
            typed_instance->m_fov = *static_cast<float*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "m_fov");
        }
        static void* get_m_fov(void* instance){ return static_cast<void*>(&(static_cast<CameraComponent*>(instance)->m_fov));}
        static bool isArray_m_fov(){ return false; }
        static const char* getFieldName_m_aspect(){ return "m_aspect";}
        static const char* getFieldTypeName_m_aspect(){ return "float";}
        static void set_m_aspect(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<CameraComponent*>(instance);
            typed_instance->m_aspect = *static_cast<float*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "m_aspect");
        }
        static void* get_m_aspect(void* instance){ return static_cast<void*>(&(static_cast<CameraComponent*>(instance)->m_aspect));}
        static bool isArray_m_aspect(){ return false; }
        static const char* getFieldName_m_zNear(){ return "m_zNear";}
        static const char* getFieldTypeName_m_zNear(){ return "float";}
        static void set_m_zNear(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<CameraComponent*>(instance);
            typed_instance->m_zNear = *static_cast<float*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "m_zNear");
        }
        static void* get_m_zNear(void* instance){ return static_cast<void*>(&(static_cast<CameraComponent*>(instance)->m_zNear));}
        static bool isArray_m_zNear(){ return false; }
        static const char* getFieldName_m_zFar(){ return "m_zFar";}
        static const char* getFieldTypeName_m_zFar(){ return "float";}
        static void set_m_zFar(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<CameraComponent*>(instance);
            typed_instance->m_zFar = *static_cast<float*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "m_zFar");
        }
        static void* get_m_zFar(void* instance){ return static_cast<void*>(&(static_cast<CameraComponent*>(instance)->m_zFar));}
        static bool isArray_m_zFar(){ return false; }
        static const char* getFieldName_isMainCamera(){ return "isMainCamera";}
        static const char* getFieldTypeName_isMainCamera(){ return "bool";}
        static void set_isMainCamera(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<CameraComponent*>(instance);
            typed_instance->isMainCamera = *static_cast<bool*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "isMainCamera");
        }
        static void* get_isMainCamera(void* instance){ return static_cast<void*>(&(static_cast<CameraComponent*>(instance)->isMainCamera));}
        static bool isArray_isMainCamera(){ return false; }
    };
}//namespace TypeFieldReflectionOparator


    void TypeWrapperRegister_CameraComponent(){
        std::unordered_map<std::string, std::string> meta_tags_m_fov;
        meta_tags_m_fov.insert({"Enable", ""});

        FieldFunctionTuple* f_field_function_tuple_m_fov=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeCameraComponentOperator::set_m_fov,
            &TypeFieldReflectionOparator::TypeCameraComponentOperator::get_m_fov,
            &TypeFieldReflectionOparator::TypeCameraComponentOperator::getClassName,
            &TypeFieldReflectionOparator::TypeCameraComponentOperator::getFieldName_m_fov,
            &TypeFieldReflectionOparator::TypeCameraComponentOperator::getFieldTypeName_m_fov,
            &TypeFieldReflectionOparator::TypeCameraComponentOperator::isArray_m_fov,
            meta_tags_m_fov
        );
        REGISTER_FIELD_TO_MAP("CameraComponent", f_field_function_tuple_m_fov);
        std::unordered_map<std::string, std::string> meta_tags_m_aspect;
        meta_tags_m_aspect.insert({"Enable", ""});

        FieldFunctionTuple* f_field_function_tuple_m_aspect=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeCameraComponentOperator::set_m_aspect,
            &TypeFieldReflectionOparator::TypeCameraComponentOperator::get_m_aspect,
            &TypeFieldReflectionOparator::TypeCameraComponentOperator::getClassName,
            &TypeFieldReflectionOparator::TypeCameraComponentOperator::getFieldName_m_aspect,
            &TypeFieldReflectionOparator::TypeCameraComponentOperator::getFieldTypeName_m_aspect,
            &TypeFieldReflectionOparator::TypeCameraComponentOperator::isArray_m_aspect,
            meta_tags_m_aspect
        );
        REGISTER_FIELD_TO_MAP("CameraComponent", f_field_function_tuple_m_aspect);
        std::unordered_map<std::string, std::string> meta_tags_m_zNear;
        meta_tags_m_zNear.insert({"Enable", ""});

        FieldFunctionTuple* f_field_function_tuple_m_zNear=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeCameraComponentOperator::set_m_zNear,
            &TypeFieldReflectionOparator::TypeCameraComponentOperator::get_m_zNear,
            &TypeFieldReflectionOparator::TypeCameraComponentOperator::getClassName,
            &TypeFieldReflectionOparator::TypeCameraComponentOperator::getFieldName_m_zNear,
            &TypeFieldReflectionOparator::TypeCameraComponentOperator::getFieldTypeName_m_zNear,
            &TypeFieldReflectionOparator::TypeCameraComponentOperator::isArray_m_zNear,
            meta_tags_m_zNear
        );
        REGISTER_FIELD_TO_MAP("CameraComponent", f_field_function_tuple_m_zNear);
        std::unordered_map<std::string, std::string> meta_tags_m_zFar;
        meta_tags_m_zFar.insert({"Enable", ""});

        FieldFunctionTuple* f_field_function_tuple_m_zFar=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeCameraComponentOperator::set_m_zFar,
            &TypeFieldReflectionOparator::TypeCameraComponentOperator::get_m_zFar,
            &TypeFieldReflectionOparator::TypeCameraComponentOperator::getClassName,
            &TypeFieldReflectionOparator::TypeCameraComponentOperator::getFieldName_m_zFar,
            &TypeFieldReflectionOparator::TypeCameraComponentOperator::getFieldTypeName_m_zFar,
            &TypeFieldReflectionOparator::TypeCameraComponentOperator::isArray_m_zFar,
            meta_tags_m_zFar
        );
        REGISTER_FIELD_TO_MAP("CameraComponent", f_field_function_tuple_m_zFar);
        std::unordered_map<std::string, std::string> meta_tags_isMainCamera;
        meta_tags_isMainCamera.insert({"Enable", ""});

        FieldFunctionTuple* f_field_function_tuple_isMainCamera=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeCameraComponentOperator::set_isMainCamera,
            &TypeFieldReflectionOparator::TypeCameraComponentOperator::get_isMainCamera,
            &TypeFieldReflectionOparator::TypeCameraComponentOperator::getClassName,
            &TypeFieldReflectionOparator::TypeCameraComponentOperator::getFieldName_isMainCamera,
            &TypeFieldReflectionOparator::TypeCameraComponentOperator::getFieldTypeName_isMainCamera,
            &TypeFieldReflectionOparator::TypeCameraComponentOperator::isArray_isMainCamera,
            meta_tags_isMainCamera
        );
        REGISTER_FIELD_TO_MAP("CameraComponent", f_field_function_tuple_isMainCamera);
        
        
        ClassFunctionTuple* f_class_function_tuple_CameraComponent=new ClassFunctionTuple(
            &TypeFieldReflectionOparator::TypeCameraComponentOperator::getCameraComponentBaseClassReflectionInstanceList,
            &TypeFieldReflectionOparator::TypeCameraComponentOperator::constructorWithJson,
            &TypeFieldReflectionOparator::TypeCameraComponentOperator::writeByName);
        REGISTER_BASE_CLASS_TO_MAP("CameraComponent", f_class_function_tuple_CameraComponent);
    }
namespace TypeWrappersRegister{
        void CameraComponent(){ TypeWrapperRegister_CameraComponent();}
}//namespace TypeWrappersRegister

}//namespace Reflection
}//namespace Lizeral


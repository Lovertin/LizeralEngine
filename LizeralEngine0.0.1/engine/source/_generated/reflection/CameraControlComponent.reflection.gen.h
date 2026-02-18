#pragma once
#include "runtime\function\ecs\components\Camera\CameraControlComponent.h"

namespace Lizeral{
    class CameraControlComponent;
namespace Reflection{
namespace TypeFieldReflectionOparator{
    class TypeCameraControlComponentOperator{
    public:
        static const char* getClassName(){ return "CameraControlComponent";}
        static void* constructorWithJson(const PJson& json_context){
            CameraControlComponent* ret_instance= new CameraControlComponent;
            PSerializer::read(json_context, *ret_instance);
            return ret_instance;
        }
        static PJson writeByName(void* instance){
            return PSerializer::write(*(CameraControlComponent*)instance);
        }
        // base class
        static int getCameraControlComponentBaseClassReflectionInstanceList(ReflectionInstance* &out_list, void* instance){
            int count = 1;
            out_list = new ReflectionInstance[count];
            for (int i=0;i<count;++i){
               out_list[i] = TypeMetaDef(Lizeral::Component,static_cast<CameraControlComponent*>(instance));
            }
            return count;
        }
        // fields
        static const char* getFieldName_move_speed(){ return "move_speed";}
        static const char* getFieldTypeName_move_speed(){ return "float";}
        static void set_move_speed(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<CameraControlComponent*>(instance);
            typed_instance->move_speed = *static_cast<float*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "move_speed");
        }
        static void* get_move_speed(void* instance){ return static_cast<void*>(&(static_cast<CameraControlComponent*>(instance)->move_speed));}
        static bool isArray_move_speed(){ return false; }
        static const char* getFieldName_m_speedMultiplier(){ return "m_speedMultiplier";}
        static const char* getFieldTypeName_m_speedMultiplier(){ return "float";}
        static void set_m_speedMultiplier(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<CameraControlComponent*>(instance);
            typed_instance->m_speedMultiplier = *static_cast<float*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "m_speedMultiplier");
        }
        static void* get_m_speedMultiplier(void* instance){ return static_cast<void*>(&(static_cast<CameraControlComponent*>(instance)->m_speedMultiplier));}
        static bool isArray_m_speedMultiplier(){ return false; }
        static const char* getFieldName_m_sensitivityX(){ return "m_sensitivityX";}
        static const char* getFieldTypeName_m_sensitivityX(){ return "float";}
        static void set_m_sensitivityX(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<CameraControlComponent*>(instance);
            typed_instance->m_sensitivityX = *static_cast<float*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "m_sensitivityX");
        }
        static void* get_m_sensitivityX(void* instance){ return static_cast<void*>(&(static_cast<CameraControlComponent*>(instance)->m_sensitivityX));}
        static bool isArray_m_sensitivityX(){ return false; }
        static const char* getFieldName_m_sensitivityY(){ return "m_sensitivityY";}
        static const char* getFieldTypeName_m_sensitivityY(){ return "float";}
        static void set_m_sensitivityY(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<CameraControlComponent*>(instance);
            typed_instance->m_sensitivityY = *static_cast<float*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "m_sensitivityY");
        }
        static void* get_m_sensitivityY(void* instance){ return static_cast<void*>(&(static_cast<CameraControlComponent*>(instance)->m_sensitivityY));}
        static bool isArray_m_sensitivityY(){ return false; }
        static const char* getFieldName_m_pitch(){ return "m_pitch";}
        static const char* getFieldTypeName_m_pitch(){ return "float";}
        static void set_m_pitch(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<CameraControlComponent*>(instance);
            typed_instance->m_pitch = *static_cast<float*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "m_pitch");
        }
        static void* get_m_pitch(void* instance){ return static_cast<void*>(&(static_cast<CameraControlComponent*>(instance)->m_pitch));}
        static bool isArray_m_pitch(){ return false; }
        static const char* getFieldName_m_yaw(){ return "m_yaw";}
        static const char* getFieldTypeName_m_yaw(){ return "float";}
        static void set_m_yaw(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<CameraControlComponent*>(instance);
            typed_instance->m_yaw = *static_cast<float*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "m_yaw");
        }
        static void* get_m_yaw(void* instance){ return static_cast<void*>(&(static_cast<CameraControlComponent*>(instance)->m_yaw));}
        static bool isArray_m_yaw(){ return false; }
    };
}//namespace TypeFieldReflectionOparator


    void TypeWrapperRegister_CameraControlComponent(){
        FieldFunctionTuple* f_field_function_tuple_move_speed=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeCameraControlComponentOperator::set_move_speed,
            &TypeFieldReflectionOparator::TypeCameraControlComponentOperator::get_move_speed,
            &TypeFieldReflectionOparator::TypeCameraControlComponentOperator::getClassName,
            &TypeFieldReflectionOparator::TypeCameraControlComponentOperator::getFieldName_move_speed,
            &TypeFieldReflectionOparator::TypeCameraControlComponentOperator::getFieldTypeName_move_speed,
            &TypeFieldReflectionOparator::TypeCameraControlComponentOperator::isArray_move_speed);
        REGISTER_FIELD_TO_MAP("CameraControlComponent", f_field_function_tuple_move_speed);
        FieldFunctionTuple* f_field_function_tuple_m_speedMultiplier=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeCameraControlComponentOperator::set_m_speedMultiplier,
            &TypeFieldReflectionOparator::TypeCameraControlComponentOperator::get_m_speedMultiplier,
            &TypeFieldReflectionOparator::TypeCameraControlComponentOperator::getClassName,
            &TypeFieldReflectionOparator::TypeCameraControlComponentOperator::getFieldName_m_speedMultiplier,
            &TypeFieldReflectionOparator::TypeCameraControlComponentOperator::getFieldTypeName_m_speedMultiplier,
            &TypeFieldReflectionOparator::TypeCameraControlComponentOperator::isArray_m_speedMultiplier);
        REGISTER_FIELD_TO_MAP("CameraControlComponent", f_field_function_tuple_m_speedMultiplier);
        FieldFunctionTuple* f_field_function_tuple_m_sensitivityX=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeCameraControlComponentOperator::set_m_sensitivityX,
            &TypeFieldReflectionOparator::TypeCameraControlComponentOperator::get_m_sensitivityX,
            &TypeFieldReflectionOparator::TypeCameraControlComponentOperator::getClassName,
            &TypeFieldReflectionOparator::TypeCameraControlComponentOperator::getFieldName_m_sensitivityX,
            &TypeFieldReflectionOparator::TypeCameraControlComponentOperator::getFieldTypeName_m_sensitivityX,
            &TypeFieldReflectionOparator::TypeCameraControlComponentOperator::isArray_m_sensitivityX);
        REGISTER_FIELD_TO_MAP("CameraControlComponent", f_field_function_tuple_m_sensitivityX);
        FieldFunctionTuple* f_field_function_tuple_m_sensitivityY=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeCameraControlComponentOperator::set_m_sensitivityY,
            &TypeFieldReflectionOparator::TypeCameraControlComponentOperator::get_m_sensitivityY,
            &TypeFieldReflectionOparator::TypeCameraControlComponentOperator::getClassName,
            &TypeFieldReflectionOparator::TypeCameraControlComponentOperator::getFieldName_m_sensitivityY,
            &TypeFieldReflectionOparator::TypeCameraControlComponentOperator::getFieldTypeName_m_sensitivityY,
            &TypeFieldReflectionOparator::TypeCameraControlComponentOperator::isArray_m_sensitivityY);
        REGISTER_FIELD_TO_MAP("CameraControlComponent", f_field_function_tuple_m_sensitivityY);
        FieldFunctionTuple* f_field_function_tuple_m_pitch=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeCameraControlComponentOperator::set_m_pitch,
            &TypeFieldReflectionOparator::TypeCameraControlComponentOperator::get_m_pitch,
            &TypeFieldReflectionOparator::TypeCameraControlComponentOperator::getClassName,
            &TypeFieldReflectionOparator::TypeCameraControlComponentOperator::getFieldName_m_pitch,
            &TypeFieldReflectionOparator::TypeCameraControlComponentOperator::getFieldTypeName_m_pitch,
            &TypeFieldReflectionOparator::TypeCameraControlComponentOperator::isArray_m_pitch);
        REGISTER_FIELD_TO_MAP("CameraControlComponent", f_field_function_tuple_m_pitch);
        FieldFunctionTuple* f_field_function_tuple_m_yaw=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeCameraControlComponentOperator::set_m_yaw,
            &TypeFieldReflectionOparator::TypeCameraControlComponentOperator::get_m_yaw,
            &TypeFieldReflectionOparator::TypeCameraControlComponentOperator::getClassName,
            &TypeFieldReflectionOparator::TypeCameraControlComponentOperator::getFieldName_m_yaw,
            &TypeFieldReflectionOparator::TypeCameraControlComponentOperator::getFieldTypeName_m_yaw,
            &TypeFieldReflectionOparator::TypeCameraControlComponentOperator::isArray_m_yaw);
        REGISTER_FIELD_TO_MAP("CameraControlComponent", f_field_function_tuple_m_yaw);
        
        
        ClassFunctionTuple* f_class_function_tuple_CameraControlComponent=new ClassFunctionTuple(
            &TypeFieldReflectionOparator::TypeCameraControlComponentOperator::getCameraControlComponentBaseClassReflectionInstanceList,
            &TypeFieldReflectionOparator::TypeCameraControlComponentOperator::constructorWithJson,
            &TypeFieldReflectionOparator::TypeCameraControlComponentOperator::writeByName);
        REGISTER_BASE_CLASS_TO_MAP("CameraControlComponent", f_class_function_tuple_CameraControlComponent);
    }
namespace TypeWrappersRegister{
        void CameraControlComponent(){ TypeWrapperRegister_CameraControlComponent();}
}//namespace TypeWrappersRegister

}//namespace Reflection
}//namespace Lizeral


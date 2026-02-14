#pragma once
#include "runtime\function\ecs\components\Collider\ColliderComponent.h"

namespace Lizeral{
    class ColliderComponent;
namespace Reflection{
namespace TypeFieldReflectionOparator{
    class TypeColliderComponentOperator{
    public:
        static const char* getClassName(){ return "ColliderComponent";}
        static void* constructorWithJson(const PJson& json_context){
            ColliderComponent* ret_instance= new ColliderComponent;
            PSerializer::read(json_context, *ret_instance);
            return ret_instance;
        }
        static PJson writeByName(void* instance){
            return PSerializer::write(*(ColliderComponent*)instance);
        }
        // base class
        static int getColliderComponentBaseClassReflectionInstanceList(ReflectionInstance* &out_list, void* instance){
            int count = 1;
            out_list = new ReflectionInstance[count];
            for (int i=0;i<count;++i){
               out_list[i] = TypeMetaDef(Lizeral::Component,static_cast<ColliderComponent*>(instance));
            }
            return count;
        }
        // fields
        static const char* getFieldName_m_type(){ return "m_type";}
        static const char* getFieldTypeName_m_type(){ return "Lizeral::ColliderType";}
        static void set_m_type(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<ColliderComponent*>(instance);
            typed_instance->m_type = *static_cast<Lizeral::ColliderType*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "m_type");
        }
        static void* get_m_type(void* instance){ return static_cast<void*>(&(static_cast<ColliderComponent*>(instance)->m_type));}
        static bool isArray_m_type(){ return false; }
        static const char* getFieldName_m_size(){ return "m_size";}
        static const char* getFieldTypeName_m_size(){ return "Lizeral::Vector3";}
        static void set_m_size(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<ColliderComponent*>(instance);
            typed_instance->m_size = *static_cast<Lizeral::Vector3*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "m_size");
        }
        static void* get_m_size(void* instance){ return static_cast<void*>(&(static_cast<ColliderComponent*>(instance)->m_size));}
        static bool isArray_m_size(){ return false; }
        static const char* getFieldName_m_radius(){ return "m_radius";}
        static const char* getFieldTypeName_m_radius(){ return "float";}
        static void set_m_radius(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<ColliderComponent*>(instance);
            typed_instance->m_radius = *static_cast<float*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "m_radius");
        }
        static void* get_m_radius(void* instance){ return static_cast<void*>(&(static_cast<ColliderComponent*>(instance)->m_radius));}
        static bool isArray_m_radius(){ return false; }
        static const char* getFieldName_m_height(){ return "m_height";}
        static const char* getFieldTypeName_m_height(){ return "float";}
        static void set_m_height(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<ColliderComponent*>(instance);
            typed_instance->m_height = *static_cast<float*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "m_height");
        }
        static void* get_m_height(void* instance){ return static_cast<void*>(&(static_cast<ColliderComponent*>(instance)->m_height));}
        static bool isArray_m_height(){ return false; }
        static const char* getFieldName_m_offset(){ return "m_offset";}
        static const char* getFieldTypeName_m_offset(){ return "Lizeral::Vector3";}
        static void set_m_offset(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<ColliderComponent*>(instance);
            typed_instance->m_offset = *static_cast<Lizeral::Vector3*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "m_offset");
        }
        static void* get_m_offset(void* instance){ return static_cast<void*>(&(static_cast<ColliderComponent*>(instance)->m_offset));}
        static bool isArray_m_offset(){ return false; }
    };
}//namespace TypeFieldReflectionOparator


    void TypeWrapperRegister_ColliderComponent(){
        FieldFunctionTuple* f_field_function_tuple_m_type=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeColliderComponentOperator::set_m_type,
            &TypeFieldReflectionOparator::TypeColliderComponentOperator::get_m_type,
            &TypeFieldReflectionOparator::TypeColliderComponentOperator::getClassName,
            &TypeFieldReflectionOparator::TypeColliderComponentOperator::getFieldName_m_type,
            &TypeFieldReflectionOparator::TypeColliderComponentOperator::getFieldTypeName_m_type,
            &TypeFieldReflectionOparator::TypeColliderComponentOperator::isArray_m_type);
        REGISTER_FIELD_TO_MAP("ColliderComponent", f_field_function_tuple_m_type);
        FieldFunctionTuple* f_field_function_tuple_m_size=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeColliderComponentOperator::set_m_size,
            &TypeFieldReflectionOparator::TypeColliderComponentOperator::get_m_size,
            &TypeFieldReflectionOparator::TypeColliderComponentOperator::getClassName,
            &TypeFieldReflectionOparator::TypeColliderComponentOperator::getFieldName_m_size,
            &TypeFieldReflectionOparator::TypeColliderComponentOperator::getFieldTypeName_m_size,
            &TypeFieldReflectionOparator::TypeColliderComponentOperator::isArray_m_size);
        REGISTER_FIELD_TO_MAP("ColliderComponent", f_field_function_tuple_m_size);
        FieldFunctionTuple* f_field_function_tuple_m_radius=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeColliderComponentOperator::set_m_radius,
            &TypeFieldReflectionOparator::TypeColliderComponentOperator::get_m_radius,
            &TypeFieldReflectionOparator::TypeColliderComponentOperator::getClassName,
            &TypeFieldReflectionOparator::TypeColliderComponentOperator::getFieldName_m_radius,
            &TypeFieldReflectionOparator::TypeColliderComponentOperator::getFieldTypeName_m_radius,
            &TypeFieldReflectionOparator::TypeColliderComponentOperator::isArray_m_radius);
        REGISTER_FIELD_TO_MAP("ColliderComponent", f_field_function_tuple_m_radius);
        FieldFunctionTuple* f_field_function_tuple_m_height=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeColliderComponentOperator::set_m_height,
            &TypeFieldReflectionOparator::TypeColliderComponentOperator::get_m_height,
            &TypeFieldReflectionOparator::TypeColliderComponentOperator::getClassName,
            &TypeFieldReflectionOparator::TypeColliderComponentOperator::getFieldName_m_height,
            &TypeFieldReflectionOparator::TypeColliderComponentOperator::getFieldTypeName_m_height,
            &TypeFieldReflectionOparator::TypeColliderComponentOperator::isArray_m_height);
        REGISTER_FIELD_TO_MAP("ColliderComponent", f_field_function_tuple_m_height);
        FieldFunctionTuple* f_field_function_tuple_m_offset=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeColliderComponentOperator::set_m_offset,
            &TypeFieldReflectionOparator::TypeColliderComponentOperator::get_m_offset,
            &TypeFieldReflectionOparator::TypeColliderComponentOperator::getClassName,
            &TypeFieldReflectionOparator::TypeColliderComponentOperator::getFieldName_m_offset,
            &TypeFieldReflectionOparator::TypeColliderComponentOperator::getFieldTypeName_m_offset,
            &TypeFieldReflectionOparator::TypeColliderComponentOperator::isArray_m_offset);
        REGISTER_FIELD_TO_MAP("ColliderComponent", f_field_function_tuple_m_offset);
        
        
        ClassFunctionTuple* f_class_function_tuple_ColliderComponent=new ClassFunctionTuple(
            &TypeFieldReflectionOparator::TypeColliderComponentOperator::getColliderComponentBaseClassReflectionInstanceList,
            &TypeFieldReflectionOparator::TypeColliderComponentOperator::constructorWithJson,
            &TypeFieldReflectionOparator::TypeColliderComponentOperator::writeByName);
        REGISTER_BASE_CLASS_TO_MAP("ColliderComponent", f_class_function_tuple_ColliderComponent);
    }
namespace TypeWrappersRegister{
        void ColliderComponent(){ TypeWrapperRegister_ColliderComponent();}
}//namespace TypeWrappersRegister

}//namespace Reflection
}//namespace Lizeral


#pragma once
#include "runtime\core\math\transform.h"

namespace Lizeral{
    class Transform;
namespace Reflection{
namespace TypeFieldReflectionOparator{
    class TypeTransformOperator{
    public:
        static const char* getClassName(){ return "Transform";}
        static void* constructorWithJson(const PJson& json_context){
            Transform* ret_instance= new Transform;
            PSerializer::read(json_context, *ret_instance);
            return ret_instance;
        }
        static PJson writeByName(void* instance){
            return PSerializer::write(*(Transform*)instance);
        }
        // base class
        static std::vector<ReflectionInstance> getTransformBaseClassReflectionInstanceList(void* instance){
            std::vector<ReflectionInstance> out_list;
            
            return out_list;
        }
        // fields
        static const char* getFieldName_m_position(){ return "m_position";}
        static const char* getFieldTypeName_m_position(){ return "Lizeral::Vector3";}
        static void set_m_position(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<Transform*>(instance);
            typed_instance->m_position = *static_cast<Lizeral::Vector3*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "m_position");
        }
        static void* get_m_position(void* instance){ return static_cast<void*>(&(static_cast<Transform*>(instance)->m_position));}
        static bool isArray_m_position(){ return false; }
        static const char* getFieldName_m_scale(){ return "m_scale";}
        static const char* getFieldTypeName_m_scale(){ return "Lizeral::Vector3";}
        static void set_m_scale(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<Transform*>(instance);
            typed_instance->m_scale = *static_cast<Lizeral::Vector3*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "m_scale");
        }
        static void* get_m_scale(void* instance){ return static_cast<void*>(&(static_cast<Transform*>(instance)->m_scale));}
        static bool isArray_m_scale(){ return false; }
        static const char* getFieldName_m_rotation(){ return "m_rotation";}
        static const char* getFieldTypeName_m_rotation(){ return "Lizeral::Quaternion";}
        static void set_m_rotation(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<Transform*>(instance);
            typed_instance->m_rotation = *static_cast<Lizeral::Quaternion*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "m_rotation");
        }
        static void* get_m_rotation(void* instance){ return static_cast<void*>(&(static_cast<Transform*>(instance)->m_rotation));}
        static bool isArray_m_rotation(){ return false; }
    };
}//namespace TypeFieldReflectionOparator


    void TypeWrapperRegister_Transform(){
        std::unordered_map<std::string, std::string> meta_tags_m_position;
        

        FieldFunctionTuple* f_field_function_tuple_m_position=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeTransformOperator::set_m_position,
            &TypeFieldReflectionOparator::TypeTransformOperator::get_m_position,
            &TypeFieldReflectionOparator::TypeTransformOperator::getClassName,
            &TypeFieldReflectionOparator::TypeTransformOperator::getFieldName_m_position,
            &TypeFieldReflectionOparator::TypeTransformOperator::getFieldTypeName_m_position,
            &TypeFieldReflectionOparator::TypeTransformOperator::isArray_m_position,
            meta_tags_m_position
        );
        REGISTER_FIELD_TO_MAP("Transform", f_field_function_tuple_m_position);
        std::unordered_map<std::string, std::string> meta_tags_m_scale;
        

        FieldFunctionTuple* f_field_function_tuple_m_scale=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeTransformOperator::set_m_scale,
            &TypeFieldReflectionOparator::TypeTransformOperator::get_m_scale,
            &TypeFieldReflectionOparator::TypeTransformOperator::getClassName,
            &TypeFieldReflectionOparator::TypeTransformOperator::getFieldName_m_scale,
            &TypeFieldReflectionOparator::TypeTransformOperator::getFieldTypeName_m_scale,
            &TypeFieldReflectionOparator::TypeTransformOperator::isArray_m_scale,
            meta_tags_m_scale
        );
        REGISTER_FIELD_TO_MAP("Transform", f_field_function_tuple_m_scale);
        std::unordered_map<std::string, std::string> meta_tags_m_rotation;
        

        FieldFunctionTuple* f_field_function_tuple_m_rotation=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeTransformOperator::set_m_rotation,
            &TypeFieldReflectionOparator::TypeTransformOperator::get_m_rotation,
            &TypeFieldReflectionOparator::TypeTransformOperator::getClassName,
            &TypeFieldReflectionOparator::TypeTransformOperator::getFieldName_m_rotation,
            &TypeFieldReflectionOparator::TypeTransformOperator::getFieldTypeName_m_rotation,
            &TypeFieldReflectionOparator::TypeTransformOperator::isArray_m_rotation,
            meta_tags_m_rotation
        );
        REGISTER_FIELD_TO_MAP("Transform", f_field_function_tuple_m_rotation);
        
        
        ClassFunctionTuple* f_class_function_tuple_Transform=new ClassFunctionTuple(
            &TypeFieldReflectionOparator::TypeTransformOperator::getTransformBaseClassReflectionInstanceList,
            &TypeFieldReflectionOparator::TypeTransformOperator::constructorWithJson,
            &TypeFieldReflectionOparator::TypeTransformOperator::writeByName);
        REGISTER_BASE_CLASS_TO_MAP("Transform", f_class_function_tuple_Transform);
    }
namespace TypeWrappersRegister{
        void Transform(){ TypeWrapperRegister_Transform();}
}//namespace TypeWrappersRegister

}//namespace Reflection
}//namespace Lizeral


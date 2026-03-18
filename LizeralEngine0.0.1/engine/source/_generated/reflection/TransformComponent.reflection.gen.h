#pragma once
#include "runtime/core/meta/serializer/serializer.h"
#include "runtime\function\ecs\components\Transform\TransformComponent.h"

namespace Lizeral{
    class TransformComponent;
namespace Reflection{
namespace TypeFieldReflectionOparator{
    class TypeTransformComponentOperator{
    public:
        static const char* getClassName(){ return "TransformComponent";}
        static void* constructorWithJson(const PJson& json_context){
            TransformComponent* ret_instance= new TransformComponent;
            PSerializer::read(json_context, *ret_instance);
            return ret_instance;
        }
        static PJson writeByName(void* instance){
            return PSerializer::write(*(TransformComponent*)instance);
        }
        // base class
        static std::vector<ReflectionInstance> getTransformComponentBaseClassReflectionInstanceList(void* instance){
            std::vector<ReflectionInstance> out_list;
            out_list.push_back(TypeMetaDef(Lizeral::Component, static_cast<TransformComponent*>(instance)));
            return out_list;
        }
        // fields
        static const char* getFieldName_m_position(){ return "m_position";}
        static const char* getFieldTypeName_m_position(){ return "Lizeral::Vector3";}
        static void set_m_position(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<TransformComponent*>(instance);
            typed_instance->m_position = *static_cast<Lizeral::Vector3*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "m_position");
        }
        static void* get_m_position(void* instance){ return static_cast<void*>(&(static_cast<TransformComponent*>(instance)->m_position));}
        static bool isArray_m_position(){ return false; }
        static const char* getFieldName_m_rotation(){ return "m_rotation";}
        static const char* getFieldTypeName_m_rotation(){ return "Lizeral::Quaternion";}
        static void set_m_rotation(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<TransformComponent*>(instance);
            typed_instance->m_rotation = *static_cast<Lizeral::Quaternion*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "m_rotation");
        }
        static void* get_m_rotation(void* instance){ return static_cast<void*>(&(static_cast<TransformComponent*>(instance)->m_rotation));}
        static bool isArray_m_rotation(){ return false; }
        static const char* getFieldName_m_scale(){ return "m_scale";}
        static const char* getFieldTypeName_m_scale(){ return "Lizeral::Vector3";}
        static void set_m_scale(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<TransformComponent*>(instance);
            typed_instance->m_scale = *static_cast<Lizeral::Vector3*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "m_scale");
        }
        static void* get_m_scale(void* instance){ return static_cast<void*>(&(static_cast<TransformComponent*>(instance)->m_scale));}
        static bool isArray_m_scale(){ return false; }
    };
}//namespace TypeFieldReflectionOparator


    void TypeWrapperRegister_TransformComponent(){
        std::unordered_map<std::string, std::string> meta_tags_m_position;
        meta_tags_m_position.insert({"Infinite", ""});
        meta_tags_m_position.insert({"UI", "Slider"});
        meta_tags_m_position.insert({"Enable", ""});

        FieldFunctionTuple* f_field_function_tuple_m_position=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeTransformComponentOperator::set_m_position,
            &TypeFieldReflectionOparator::TypeTransformComponentOperator::get_m_position,
            &TypeFieldReflectionOparator::TypeTransformComponentOperator::getClassName,
            &TypeFieldReflectionOparator::TypeTransformComponentOperator::getFieldName_m_position,
            &TypeFieldReflectionOparator::TypeTransformComponentOperator::getFieldTypeName_m_position,
            &TypeFieldReflectionOparator::TypeTransformComponentOperator::isArray_m_position,
            meta_tags_m_position
        );
        REGISTER_FIELD_TO_MAP("TransformComponent", f_field_function_tuple_m_position);
        std::unordered_map<std::string, std::string> meta_tags_m_rotation;
        meta_tags_m_rotation.insert({"Range", "-360~360"});
        meta_tags_m_rotation.insert({"UI", "Slider"});
        meta_tags_m_rotation.insert({"Enable", ""});

        FieldFunctionTuple* f_field_function_tuple_m_rotation=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeTransformComponentOperator::set_m_rotation,
            &TypeFieldReflectionOparator::TypeTransformComponentOperator::get_m_rotation,
            &TypeFieldReflectionOparator::TypeTransformComponentOperator::getClassName,
            &TypeFieldReflectionOparator::TypeTransformComponentOperator::getFieldName_m_rotation,
            &TypeFieldReflectionOparator::TypeTransformComponentOperator::getFieldTypeName_m_rotation,
            &TypeFieldReflectionOparator::TypeTransformComponentOperator::isArray_m_rotation,
            meta_tags_m_rotation
        );
        REGISTER_FIELD_TO_MAP("TransformComponent", f_field_function_tuple_m_rotation);
        std::unordered_map<std::string, std::string> meta_tags_m_scale;
        meta_tags_m_scale.insert({"Enable", ""});

        FieldFunctionTuple* f_field_function_tuple_m_scale=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeTransformComponentOperator::set_m_scale,
            &TypeFieldReflectionOparator::TypeTransformComponentOperator::get_m_scale,
            &TypeFieldReflectionOparator::TypeTransformComponentOperator::getClassName,
            &TypeFieldReflectionOparator::TypeTransformComponentOperator::getFieldName_m_scale,
            &TypeFieldReflectionOparator::TypeTransformComponentOperator::getFieldTypeName_m_scale,
            &TypeFieldReflectionOparator::TypeTransformComponentOperator::isArray_m_scale,
            meta_tags_m_scale
        );
        REGISTER_FIELD_TO_MAP("TransformComponent", f_field_function_tuple_m_scale);
        
        
        ClassFunctionTuple* f_class_function_tuple_TransformComponent=new ClassFunctionTuple(
            &TypeFieldReflectionOparator::TypeTransformComponentOperator::getTransformComponentBaseClassReflectionInstanceList,
            &TypeFieldReflectionOparator::TypeTransformComponentOperator::constructorWithJson,
            &TypeFieldReflectionOparator::TypeTransformComponentOperator::writeByName);
        REGISTER_BASE_CLASS_TO_MAP("TransformComponent", f_class_function_tuple_TransformComponent);
    }
namespace TypeWrappersRegister{
        void TransformComponent(){ TypeWrapperRegister_TransformComponent();}
}//namespace TypeWrappersRegister

}//namespace Reflection
}//namespace Lizeral


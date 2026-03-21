#pragma once
#include "runtime/core/meta/serializer/serializer.h"
#include "runtime\function\ecs\components\Light\DirectionalLightComponent.h"

namespace Lizeral{
    class DirectionLightComponent;
namespace Reflection{
namespace TypeFieldReflectionOparator{
    class TypeDirectionLightComponentOperator{
    public:
        static const char* getClassName(){ return "DirectionLightComponent";}
        static void* constructorWithJson(const PJson& json_context){
            DirectionLightComponent* ret_instance= new DirectionLightComponent;
            PSerializer::read(json_context, *ret_instance);
            return ret_instance;
        }
        static PJson writeByName(void* instance){
            return PSerializer::write(*(DirectionLightComponent*)instance);
        }
        // base class
        static std::vector<ReflectionInstance> getDirectionLightComponentBaseClassReflectionInstanceList(void* instance){
            std::vector<ReflectionInstance> out_list;
            out_list.push_back(TypeMetaDef(Lizeral::Component, static_cast<DirectionLightComponent*>(instance)));
            return out_list;
        }
        // fields
        static const char* getFieldName_m_color(){ return "m_color";}
        static const char* getFieldTypeName_m_color(){ return "Lizeral::Vector3";}
        static void set_m_color(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<DirectionLightComponent*>(instance);
            typed_instance->m_color = *static_cast<Lizeral::Vector3*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "m_color");
        }
        static void* get_m_color(void* instance){ return static_cast<void*>(&(static_cast<DirectionLightComponent*>(instance)->m_color));}
        static bool isArray_m_color(){ return false; }
        static const char* getFieldName_m_intensity(){ return "m_intensity";}
        static const char* getFieldTypeName_m_intensity(){ return "float";}
        static void set_m_intensity(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<DirectionLightComponent*>(instance);
            typed_instance->m_intensity = *static_cast<float*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "m_intensity");
        }
        static void* get_m_intensity(void* instance){ return static_cast<void*>(&(static_cast<DirectionLightComponent*>(instance)->m_intensity));}
        static bool isArray_m_intensity(){ return false; }
        static const char* getFieldName_m_isGlobal(){ return "m_isGlobal";}
        static const char* getFieldTypeName_m_isGlobal(){ return "bool";}
        static void set_m_isGlobal(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<DirectionLightComponent*>(instance);
            typed_instance->m_isGlobal = *static_cast<bool*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "m_isGlobal");
        }
        static void* get_m_isGlobal(void* instance){ return static_cast<void*>(&(static_cast<DirectionLightComponent*>(instance)->m_isGlobal));}
        static bool isArray_m_isGlobal(){ return false; }
    };
}//namespace TypeFieldReflectionOparator


    void TypeWrapperRegister_DirectionLightComponent(){
        std::unordered_map<std::string, std::string> meta_tags_m_color;
        meta_tags_m_color.insert({"Enable", ""});

        FieldFunctionTuple* f_field_function_tuple_m_color=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeDirectionLightComponentOperator::set_m_color,
            &TypeFieldReflectionOparator::TypeDirectionLightComponentOperator::get_m_color,
            &TypeFieldReflectionOparator::TypeDirectionLightComponentOperator::getClassName,
            &TypeFieldReflectionOparator::TypeDirectionLightComponentOperator::getFieldName_m_color,
            &TypeFieldReflectionOparator::TypeDirectionLightComponentOperator::getFieldTypeName_m_color,
            &TypeFieldReflectionOparator::TypeDirectionLightComponentOperator::isArray_m_color,
            meta_tags_m_color
        );
        REGISTER_FIELD_TO_MAP("DirectionLightComponent", f_field_function_tuple_m_color);
        std::unordered_map<std::string, std::string> meta_tags_m_intensity;
        meta_tags_m_intensity.insert({"Enable", ""});

        FieldFunctionTuple* f_field_function_tuple_m_intensity=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeDirectionLightComponentOperator::set_m_intensity,
            &TypeFieldReflectionOparator::TypeDirectionLightComponentOperator::get_m_intensity,
            &TypeFieldReflectionOparator::TypeDirectionLightComponentOperator::getClassName,
            &TypeFieldReflectionOparator::TypeDirectionLightComponentOperator::getFieldName_m_intensity,
            &TypeFieldReflectionOparator::TypeDirectionLightComponentOperator::getFieldTypeName_m_intensity,
            &TypeFieldReflectionOparator::TypeDirectionLightComponentOperator::isArray_m_intensity,
            meta_tags_m_intensity
        );
        REGISTER_FIELD_TO_MAP("DirectionLightComponent", f_field_function_tuple_m_intensity);
        std::unordered_map<std::string, std::string> meta_tags_m_isGlobal;
        meta_tags_m_isGlobal.insert({"Enable", ""});

        FieldFunctionTuple* f_field_function_tuple_m_isGlobal=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeDirectionLightComponentOperator::set_m_isGlobal,
            &TypeFieldReflectionOparator::TypeDirectionLightComponentOperator::get_m_isGlobal,
            &TypeFieldReflectionOparator::TypeDirectionLightComponentOperator::getClassName,
            &TypeFieldReflectionOparator::TypeDirectionLightComponentOperator::getFieldName_m_isGlobal,
            &TypeFieldReflectionOparator::TypeDirectionLightComponentOperator::getFieldTypeName_m_isGlobal,
            &TypeFieldReflectionOparator::TypeDirectionLightComponentOperator::isArray_m_isGlobal,
            meta_tags_m_isGlobal
        );
        REGISTER_FIELD_TO_MAP("DirectionLightComponent", f_field_function_tuple_m_isGlobal);
        
        
        ClassFunctionTuple* f_class_function_tuple_DirectionLightComponent=new ClassFunctionTuple(
            &TypeFieldReflectionOparator::TypeDirectionLightComponentOperator::getDirectionLightComponentBaseClassReflectionInstanceList,
            &TypeFieldReflectionOparator::TypeDirectionLightComponentOperator::constructorWithJson,
            &TypeFieldReflectionOparator::TypeDirectionLightComponentOperator::writeByName);
        REGISTER_BASE_CLASS_TO_MAP("DirectionLightComponent", f_class_function_tuple_DirectionLightComponent);
    }
namespace TypeWrappersRegister{
        void DirectionLightComponent(){ TypeWrapperRegister_DirectionLightComponent();}
}//namespace TypeWrappersRegister

}//namespace Reflection
}//namespace Lizeral


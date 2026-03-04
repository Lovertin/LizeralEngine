#pragma once
#include "runtime\function\ecs\components\Name\NameComponent.h"

namespace Lizeral{
    class NameComponent;
namespace Reflection{
namespace TypeFieldReflectionOparator{
    class TypeNameComponentOperator{
    public:
        static const char* getClassName(){ return "NameComponent";}
        static void* constructorWithJson(const PJson& json_context){
            NameComponent* ret_instance= new NameComponent;
            PSerializer::read(json_context, *ret_instance);
            return ret_instance;
        }
        static PJson writeByName(void* instance){
            return PSerializer::write(*(NameComponent*)instance);
        }
        // base class
        static std::vector<ReflectionInstance> getNameComponentBaseClassReflectionInstanceList(void* instance){
            std::vector<ReflectionInstance> out_list;
            out_list.push_back(TypeMetaDef(Lizeral::Component, static_cast<NameComponent*>(instance)));
            return out_list;
        }
        // fields
        static const char* getFieldName_m_name(){ return "m_name";}
        static const char* getFieldTypeName_m_name(){ return "std::string";}
        static void set_m_name(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<NameComponent*>(instance);
            typed_instance->m_name = *static_cast<std::string*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "m_name");
        }
        static void* get_m_name(void* instance){ return static_cast<void*>(&(static_cast<NameComponent*>(instance)->m_name));}
        static bool isArray_m_name(){ return false; }
    };
}//namespace TypeFieldReflectionOparator


    void TypeWrapperRegister_NameComponent(){
        std::unordered_map<std::string, std::string> meta_tags_m_name;
        meta_tags_m_name.insert({"UI", "Headline"});
        meta_tags_m_name.insert({"Enable", ""});

        FieldFunctionTuple* f_field_function_tuple_m_name=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeNameComponentOperator::set_m_name,
            &TypeFieldReflectionOparator::TypeNameComponentOperator::get_m_name,
            &TypeFieldReflectionOparator::TypeNameComponentOperator::getClassName,
            &TypeFieldReflectionOparator::TypeNameComponentOperator::getFieldName_m_name,
            &TypeFieldReflectionOparator::TypeNameComponentOperator::getFieldTypeName_m_name,
            &TypeFieldReflectionOparator::TypeNameComponentOperator::isArray_m_name,
            meta_tags_m_name
        );
        REGISTER_FIELD_TO_MAP("NameComponent", f_field_function_tuple_m_name);
        
        
        ClassFunctionTuple* f_class_function_tuple_NameComponent=new ClassFunctionTuple(
            &TypeFieldReflectionOparator::TypeNameComponentOperator::getNameComponentBaseClassReflectionInstanceList,
            &TypeFieldReflectionOparator::TypeNameComponentOperator::constructorWithJson,
            &TypeFieldReflectionOparator::TypeNameComponentOperator::writeByName);
        REGISTER_BASE_CLASS_TO_MAP("NameComponent", f_class_function_tuple_NameComponent);
    }
namespace TypeWrappersRegister{
        void NameComponent(){ TypeWrapperRegister_NameComponent();}
}//namespace TypeWrappersRegister

}//namespace Reflection
}//namespace Lizeral


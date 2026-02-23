#pragma once
#include "runtime\resource\resource.h"

namespace Lizeral{
    class Resource;
namespace Reflection{
namespace TypeFieldReflectionOparator{
    class TypeResourceOperator{
    public:
        static const char* getClassName(){ return "Resource";}
        static void* constructorWithJson(const PJson& json_context){
            Resource* ret_instance= new Resource;
            PSerializer::read(json_context, *ret_instance);
            return ret_instance;
        }
        static PJson writeByName(void* instance){
            return PSerializer::write(*(Resource*)instance);
        }
        // base class
        static int getResourceBaseClassReflectionInstanceList(ReflectionInstance* &out_list, void* instance){
            int count = 0;
            
            return count;
        }
        // fields
        static const char* getFieldName_m_path(){ return "m_path";}
        static const char* getFieldTypeName_m_path(){ return "std::string";}
        static void set_m_path(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<Resource*>(instance);
            typed_instance->m_path = *static_cast<std::string*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "m_path");
        }
        static void* get_m_path(void* instance){ return static_cast<void*>(&(static_cast<Resource*>(instance)->m_path));}
        static bool isArray_m_path(){ return false; }
    };
}//namespace TypeFieldReflectionOparator


    void TypeWrapperRegister_Resource(){
        FieldFunctionTuple* f_field_function_tuple_m_path=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeResourceOperator::set_m_path,
            &TypeFieldReflectionOparator::TypeResourceOperator::get_m_path,
            &TypeFieldReflectionOparator::TypeResourceOperator::getClassName,
            &TypeFieldReflectionOparator::TypeResourceOperator::getFieldName_m_path,
            &TypeFieldReflectionOparator::TypeResourceOperator::getFieldTypeName_m_path,
            &TypeFieldReflectionOparator::TypeResourceOperator::isArray_m_path);
        REGISTER_FIELD_TO_MAP("Resource", f_field_function_tuple_m_path);
        
        
        ClassFunctionTuple* f_class_function_tuple_Resource=new ClassFunctionTuple(
            &TypeFieldReflectionOparator::TypeResourceOperator::getResourceBaseClassReflectionInstanceList,
            &TypeFieldReflectionOparator::TypeResourceOperator::constructorWithJson,
            &TypeFieldReflectionOparator::TypeResourceOperator::writeByName);
        REGISTER_BASE_CLASS_TO_MAP("Resource", f_class_function_tuple_Resource);
    }
namespace TypeWrappersRegister{
        void Resource(){ TypeWrapperRegister_Resource();}
}//namespace TypeWrappersRegister

}//namespace Reflection
}//namespace Lizeral


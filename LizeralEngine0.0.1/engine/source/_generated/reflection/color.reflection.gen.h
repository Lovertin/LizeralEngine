#pragma once
#include "runtime/core/meta/serializer/serializer.h"
#include "runtime\core\color\color.h"

namespace Lizeral{
    class Color;
namespace Reflection{
namespace TypeFieldReflectionOparator{
    class TypeColorOperator{
    public:
        static const char* getClassName(){ return "Color";}
        static void* constructorWithJson(const PJson& json_context){
            Color* ret_instance= new Color;
            PSerializer::read(json_context, *ret_instance);
            return ret_instance;
        }
        static PJson writeByName(void* instance){
            return PSerializer::write(*(Color*)instance);
        }
        // base class
        static std::vector<ReflectionInstance> getColorBaseClassReflectionInstanceList(void* instance){
            std::vector<ReflectionInstance> out_list;
            
            return out_list;
        }
        // fields
        static const char* getFieldName_r(){ return "r";}
        static const char* getFieldTypeName_r(){ return "float";}
        static void set_r(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<Color*>(instance);
            typed_instance->r = *static_cast<float*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "r");
        }
        static void* get_r(void* instance){ return static_cast<void*>(&(static_cast<Color*>(instance)->r));}
        static bool isArray_r(){ return false; }
        static const char* getFieldName_g(){ return "g";}
        static const char* getFieldTypeName_g(){ return "float";}
        static void set_g(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<Color*>(instance);
            typed_instance->g = *static_cast<float*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "g");
        }
        static void* get_g(void* instance){ return static_cast<void*>(&(static_cast<Color*>(instance)->g));}
        static bool isArray_g(){ return false; }
        static const char* getFieldName_b(){ return "b";}
        static const char* getFieldTypeName_b(){ return "float";}
        static void set_b(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<Color*>(instance);
            typed_instance->b = *static_cast<float*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "b");
        }
        static void* get_b(void* instance){ return static_cast<void*>(&(static_cast<Color*>(instance)->b));}
        static bool isArray_b(){ return false; }
    };
}//namespace TypeFieldReflectionOparator


    void TypeWrapperRegister_Color(){
        std::unordered_map<std::string, std::string> meta_tags_r;
        

        FieldFunctionTuple* f_field_function_tuple_r=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeColorOperator::set_r,
            &TypeFieldReflectionOparator::TypeColorOperator::get_r,
            &TypeFieldReflectionOparator::TypeColorOperator::getClassName,
            &TypeFieldReflectionOparator::TypeColorOperator::getFieldName_r,
            &TypeFieldReflectionOparator::TypeColorOperator::getFieldTypeName_r,
            &TypeFieldReflectionOparator::TypeColorOperator::isArray_r,
            meta_tags_r
        );
        REGISTER_FIELD_TO_MAP("Color", f_field_function_tuple_r);
        std::unordered_map<std::string, std::string> meta_tags_g;
        

        FieldFunctionTuple* f_field_function_tuple_g=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeColorOperator::set_g,
            &TypeFieldReflectionOparator::TypeColorOperator::get_g,
            &TypeFieldReflectionOparator::TypeColorOperator::getClassName,
            &TypeFieldReflectionOparator::TypeColorOperator::getFieldName_g,
            &TypeFieldReflectionOparator::TypeColorOperator::getFieldTypeName_g,
            &TypeFieldReflectionOparator::TypeColorOperator::isArray_g,
            meta_tags_g
        );
        REGISTER_FIELD_TO_MAP("Color", f_field_function_tuple_g);
        std::unordered_map<std::string, std::string> meta_tags_b;
        

        FieldFunctionTuple* f_field_function_tuple_b=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeColorOperator::set_b,
            &TypeFieldReflectionOparator::TypeColorOperator::get_b,
            &TypeFieldReflectionOparator::TypeColorOperator::getClassName,
            &TypeFieldReflectionOparator::TypeColorOperator::getFieldName_b,
            &TypeFieldReflectionOparator::TypeColorOperator::getFieldTypeName_b,
            &TypeFieldReflectionOparator::TypeColorOperator::isArray_b,
            meta_tags_b
        );
        REGISTER_FIELD_TO_MAP("Color", f_field_function_tuple_b);
        
        
        ClassFunctionTuple* f_class_function_tuple_Color=new ClassFunctionTuple(
            &TypeFieldReflectionOparator::TypeColorOperator::getColorBaseClassReflectionInstanceList,
            &TypeFieldReflectionOparator::TypeColorOperator::constructorWithJson,
            &TypeFieldReflectionOparator::TypeColorOperator::writeByName);
        REGISTER_BASE_CLASS_TO_MAP("Color", f_class_function_tuple_Color);
    }
namespace TypeWrappersRegister{
        void Color(){ TypeWrapperRegister_Color();}
}//namespace TypeWrappersRegister

}//namespace Reflection
}//namespace Lizeral


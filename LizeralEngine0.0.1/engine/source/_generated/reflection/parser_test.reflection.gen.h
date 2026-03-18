#pragma once
#include "runtime/core/meta/serializer/serializer.h"
#include "runtime\core\test\parser_test.h"

namespace Lizeral{
    class Test;
namespace Reflection{
namespace TypeFieldReflectionOparator{
    class TypeTestOperator{
    public:
        static const char* getClassName(){ return "Test";}
        static void* constructorWithJson(const PJson& json_context){
            Test* ret_instance= new Test;
            PSerializer::read(json_context, *ret_instance);
            return ret_instance;
        }
        static PJson writeByName(void* instance){
            return PSerializer::write(*(Test*)instance);
        }
        // base class
        static std::vector<ReflectionInstance> getTestBaseClassReflectionInstanceList(void* instance){
            std::vector<ReflectionInstance> out_list;
            
            return out_list;
        }
        // fields
        static const char* getFieldName_x(){ return "x";}
        static const char* getFieldTypeName_x(){ return "float";}
        static void set_x(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<Test*>(instance);
            typed_instance->x = *static_cast<float*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "x");
        }
        static void* get_x(void* instance){ return static_cast<void*>(&(static_cast<Test*>(instance)->x));}
        static bool isArray_x(){ return false; }
        static const char* getFieldName_y(){ return "y";}
        static const char* getFieldTypeName_y(){ return "float";}
        static void set_y(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<Test*>(instance);
            typed_instance->y = *static_cast<float*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "y");
        }
        static void* get_y(void* instance){ return static_cast<void*>(&(static_cast<Test*>(instance)->y));}
        static bool isArray_y(){ return false; }
        static const char* getFieldName_z(){ return "z";}
        static const char* getFieldTypeName_z(){ return "int";}
        static void set_z(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<Test*>(instance);
            typed_instance->z = *static_cast<int*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "z");
        }
        static void* get_z(void* instance){ return static_cast<void*>(&(static_cast<Test*>(instance)->z));}
        static bool isArray_z(){ return false; }
        static const char* getFieldName_s(){ return "s";}
        static const char* getFieldTypeName_s(){ return "std::string";}
        static void set_s(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<Test*>(instance);
            typed_instance->s = *static_cast<std::string*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "s");
        }
        static void* get_s(void* instance){ return static_cast<void*>(&(static_cast<Test*>(instance)->s));}
        static bool isArray_s(){ return false; }
    };
}//namespace TypeFieldReflectionOparator


    void TypeWrapperRegister_Test(){
        std::unordered_map<std::string, std::string> meta_tags_x;
        

        FieldFunctionTuple* f_field_function_tuple_x=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeTestOperator::set_x,
            &TypeFieldReflectionOparator::TypeTestOperator::get_x,
            &TypeFieldReflectionOparator::TypeTestOperator::getClassName,
            &TypeFieldReflectionOparator::TypeTestOperator::getFieldName_x,
            &TypeFieldReflectionOparator::TypeTestOperator::getFieldTypeName_x,
            &TypeFieldReflectionOparator::TypeTestOperator::isArray_x,
            meta_tags_x
        );
        REGISTER_FIELD_TO_MAP("Test", f_field_function_tuple_x);
        std::unordered_map<std::string, std::string> meta_tags_y;
        

        FieldFunctionTuple* f_field_function_tuple_y=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeTestOperator::set_y,
            &TypeFieldReflectionOparator::TypeTestOperator::get_y,
            &TypeFieldReflectionOparator::TypeTestOperator::getClassName,
            &TypeFieldReflectionOparator::TypeTestOperator::getFieldName_y,
            &TypeFieldReflectionOparator::TypeTestOperator::getFieldTypeName_y,
            &TypeFieldReflectionOparator::TypeTestOperator::isArray_y,
            meta_tags_y
        );
        REGISTER_FIELD_TO_MAP("Test", f_field_function_tuple_y);
        std::unordered_map<std::string, std::string> meta_tags_z;
        

        FieldFunctionTuple* f_field_function_tuple_z=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeTestOperator::set_z,
            &TypeFieldReflectionOparator::TypeTestOperator::get_z,
            &TypeFieldReflectionOparator::TypeTestOperator::getClassName,
            &TypeFieldReflectionOparator::TypeTestOperator::getFieldName_z,
            &TypeFieldReflectionOparator::TypeTestOperator::getFieldTypeName_z,
            &TypeFieldReflectionOparator::TypeTestOperator::isArray_z,
            meta_tags_z
        );
        REGISTER_FIELD_TO_MAP("Test", f_field_function_tuple_z);
        std::unordered_map<std::string, std::string> meta_tags_s;
        

        FieldFunctionTuple* f_field_function_tuple_s=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeTestOperator::set_s,
            &TypeFieldReflectionOparator::TypeTestOperator::get_s,
            &TypeFieldReflectionOparator::TypeTestOperator::getClassName,
            &TypeFieldReflectionOparator::TypeTestOperator::getFieldName_s,
            &TypeFieldReflectionOparator::TypeTestOperator::getFieldTypeName_s,
            &TypeFieldReflectionOparator::TypeTestOperator::isArray_s,
            meta_tags_s
        );
        REGISTER_FIELD_TO_MAP("Test", f_field_function_tuple_s);
        
        
        ClassFunctionTuple* f_class_function_tuple_Test=new ClassFunctionTuple(
            &TypeFieldReflectionOparator::TypeTestOperator::getTestBaseClassReflectionInstanceList,
            &TypeFieldReflectionOparator::TypeTestOperator::constructorWithJson,
            &TypeFieldReflectionOparator::TypeTestOperator::writeByName);
        REGISTER_BASE_CLASS_TO_MAP("Test", f_class_function_tuple_Test);
    }
namespace TypeWrappersRegister{
        void Test(){ TypeWrapperRegister_Test();}
}//namespace TypeWrappersRegister

}//namespace Reflection
}//namespace Lizeral


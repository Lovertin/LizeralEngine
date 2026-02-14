#pragma once
#include "runtime\core\math\vector2.h"

namespace Lizeral{
    class Vector2;
namespace Reflection{
namespace TypeFieldReflectionOparator{
    class TypeVector2Operator{
    public:
        static const char* getClassName(){ return "Vector2";}
        static void* constructorWithJson(const PJson& json_context){
            Vector2* ret_instance= new Vector2;
            PSerializer::read(json_context, *ret_instance);
            return ret_instance;
        }
        static PJson writeByName(void* instance){
            return PSerializer::write(*(Vector2*)instance);
        }
        // base class
        static int getVector2BaseClassReflectionInstanceList(ReflectionInstance* &out_list, void* instance){
            int count = 0;
            
            return count;
        }
        // fields
        static const char* getFieldName_x(){ return "x";}
        static const char* getFieldTypeName_x(){ return "float";}
        static void set_x(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<Vector2*>(instance);
            typed_instance->x = *static_cast<float*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "x");
        }
        static void* get_x(void* instance){ return static_cast<void*>(&(static_cast<Vector2*>(instance)->x));}
        static bool isArray_x(){ return false; }
        static const char* getFieldName_y(){ return "y";}
        static const char* getFieldTypeName_y(){ return "float";}
        static void set_y(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<Vector2*>(instance);
            typed_instance->y = *static_cast<float*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "y");
        }
        static void* get_y(void* instance){ return static_cast<void*>(&(static_cast<Vector2*>(instance)->y));}
        static bool isArray_y(){ return false; }
    };
}//namespace TypeFieldReflectionOparator


    void TypeWrapperRegister_Vector2(){
        FieldFunctionTuple* f_field_function_tuple_x=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeVector2Operator::set_x,
            &TypeFieldReflectionOparator::TypeVector2Operator::get_x,
            &TypeFieldReflectionOparator::TypeVector2Operator::getClassName,
            &TypeFieldReflectionOparator::TypeVector2Operator::getFieldName_x,
            &TypeFieldReflectionOparator::TypeVector2Operator::getFieldTypeName_x,
            &TypeFieldReflectionOparator::TypeVector2Operator::isArray_x);
        REGISTER_FIELD_TO_MAP("Vector2", f_field_function_tuple_x);
        FieldFunctionTuple* f_field_function_tuple_y=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeVector2Operator::set_y,
            &TypeFieldReflectionOparator::TypeVector2Operator::get_y,
            &TypeFieldReflectionOparator::TypeVector2Operator::getClassName,
            &TypeFieldReflectionOparator::TypeVector2Operator::getFieldName_y,
            &TypeFieldReflectionOparator::TypeVector2Operator::getFieldTypeName_y,
            &TypeFieldReflectionOparator::TypeVector2Operator::isArray_y);
        REGISTER_FIELD_TO_MAP("Vector2", f_field_function_tuple_y);
        
        
        ClassFunctionTuple* f_class_function_tuple_Vector2=new ClassFunctionTuple(
            &TypeFieldReflectionOparator::TypeVector2Operator::getVector2BaseClassReflectionInstanceList,
            &TypeFieldReflectionOparator::TypeVector2Operator::constructorWithJson,
            &TypeFieldReflectionOparator::TypeVector2Operator::writeByName);
        REGISTER_BASE_CLASS_TO_MAP("Vector2", f_class_function_tuple_Vector2);
    }
namespace TypeWrappersRegister{
        void Vector2(){ TypeWrapperRegister_Vector2();}
}//namespace TypeWrappersRegister

}//namespace Reflection
}//namespace Lizeral


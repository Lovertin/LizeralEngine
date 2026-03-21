#pragma once
#include "runtime/core/meta/serializer/serializer.h"
#include "runtime\core\math\vector4.h"

namespace Lizeral{
    class Vector4;
namespace Reflection{
namespace TypeFieldReflectionOparator{
    class TypeVector4Operator{
    public:
        static const char* getClassName(){ return "Vector4";}
        static void* constructorWithJson(const PJson& json_context){
            Vector4* ret_instance= new Vector4;
            PSerializer::read(json_context, *ret_instance);
            return ret_instance;
        }
        static PJson writeByName(void* instance){
            return PSerializer::write(*(Vector4*)instance);
        }
        // base class
        static std::vector<ReflectionInstance> getVector4BaseClassReflectionInstanceList(void* instance){
            std::vector<ReflectionInstance> out_list;
            
            return out_list;
        }
        // fields
        static const char* getFieldName_x(){ return "x";}
        static const char* getFieldTypeName_x(){ return "float";}
        static void set_x(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<Vector4*>(instance);
            typed_instance->x = *static_cast<float*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "x");
        }
        static void* get_x(void* instance){ return static_cast<void*>(&(static_cast<Vector4*>(instance)->x));}
        static bool isArray_x(){ return false; }
        static const char* getFieldName_y(){ return "y";}
        static const char* getFieldTypeName_y(){ return "float";}
        static void set_y(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<Vector4*>(instance);
            typed_instance->y = *static_cast<float*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "y");
        }
        static void* get_y(void* instance){ return static_cast<void*>(&(static_cast<Vector4*>(instance)->y));}
        static bool isArray_y(){ return false; }
        static const char* getFieldName_z(){ return "z";}
        static const char* getFieldTypeName_z(){ return "float";}
        static void set_z(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<Vector4*>(instance);
            typed_instance->z = *static_cast<float*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "z");
        }
        static void* get_z(void* instance){ return static_cast<void*>(&(static_cast<Vector4*>(instance)->z));}
        static bool isArray_z(){ return false; }
        static const char* getFieldName_w(){ return "w";}
        static const char* getFieldTypeName_w(){ return "float";}
        static void set_w(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<Vector4*>(instance);
            typed_instance->w = *static_cast<float*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "w");
        }
        static void* get_w(void* instance){ return static_cast<void*>(&(static_cast<Vector4*>(instance)->w));}
        static bool isArray_w(){ return false; }
    };
}//namespace TypeFieldReflectionOparator


    void TypeWrapperRegister_Vector4(){
        std::unordered_map<std::string, std::string> meta_tags_x;
        

        FieldFunctionTuple* f_field_function_tuple_x=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeVector4Operator::set_x,
            &TypeFieldReflectionOparator::TypeVector4Operator::get_x,
            &TypeFieldReflectionOparator::TypeVector4Operator::getClassName,
            &TypeFieldReflectionOparator::TypeVector4Operator::getFieldName_x,
            &TypeFieldReflectionOparator::TypeVector4Operator::getFieldTypeName_x,
            &TypeFieldReflectionOparator::TypeVector4Operator::isArray_x,
            meta_tags_x
        );
        REGISTER_FIELD_TO_MAP("Vector4", f_field_function_tuple_x);
        std::unordered_map<std::string, std::string> meta_tags_y;
        

        FieldFunctionTuple* f_field_function_tuple_y=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeVector4Operator::set_y,
            &TypeFieldReflectionOparator::TypeVector4Operator::get_y,
            &TypeFieldReflectionOparator::TypeVector4Operator::getClassName,
            &TypeFieldReflectionOparator::TypeVector4Operator::getFieldName_y,
            &TypeFieldReflectionOparator::TypeVector4Operator::getFieldTypeName_y,
            &TypeFieldReflectionOparator::TypeVector4Operator::isArray_y,
            meta_tags_y
        );
        REGISTER_FIELD_TO_MAP("Vector4", f_field_function_tuple_y);
        std::unordered_map<std::string, std::string> meta_tags_z;
        

        FieldFunctionTuple* f_field_function_tuple_z=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeVector4Operator::set_z,
            &TypeFieldReflectionOparator::TypeVector4Operator::get_z,
            &TypeFieldReflectionOparator::TypeVector4Operator::getClassName,
            &TypeFieldReflectionOparator::TypeVector4Operator::getFieldName_z,
            &TypeFieldReflectionOparator::TypeVector4Operator::getFieldTypeName_z,
            &TypeFieldReflectionOparator::TypeVector4Operator::isArray_z,
            meta_tags_z
        );
        REGISTER_FIELD_TO_MAP("Vector4", f_field_function_tuple_z);
        std::unordered_map<std::string, std::string> meta_tags_w;
        

        FieldFunctionTuple* f_field_function_tuple_w=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeVector4Operator::set_w,
            &TypeFieldReflectionOparator::TypeVector4Operator::get_w,
            &TypeFieldReflectionOparator::TypeVector4Operator::getClassName,
            &TypeFieldReflectionOparator::TypeVector4Operator::getFieldName_w,
            &TypeFieldReflectionOparator::TypeVector4Operator::getFieldTypeName_w,
            &TypeFieldReflectionOparator::TypeVector4Operator::isArray_w,
            meta_tags_w
        );
        REGISTER_FIELD_TO_MAP("Vector4", f_field_function_tuple_w);
        
        
        ClassFunctionTuple* f_class_function_tuple_Vector4=new ClassFunctionTuple(
            &TypeFieldReflectionOparator::TypeVector4Operator::getVector4BaseClassReflectionInstanceList,
            &TypeFieldReflectionOparator::TypeVector4Operator::constructorWithJson,
            &TypeFieldReflectionOparator::TypeVector4Operator::writeByName);
        REGISTER_BASE_CLASS_TO_MAP("Vector4", f_class_function_tuple_Vector4);
    }
namespace TypeWrappersRegister{
        void Vector4(){ TypeWrapperRegister_Vector4();}
}//namespace TypeWrappersRegister

}//namespace Reflection
}//namespace Lizeral


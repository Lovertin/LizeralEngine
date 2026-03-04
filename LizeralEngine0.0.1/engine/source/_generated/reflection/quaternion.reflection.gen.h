#pragma once
#include "runtime\core\math\quaternion.h"

namespace Lizeral{
    class Quaternion;
namespace Reflection{
namespace TypeFieldReflectionOparator{
    class TypeQuaternionOperator{
    public:
        static const char* getClassName(){ return "Quaternion";}
        static void* constructorWithJson(const PJson& json_context){
            Quaternion* ret_instance= new Quaternion;
            PSerializer::read(json_context, *ret_instance);
            return ret_instance;
        }
        static PJson writeByName(void* instance){
            return PSerializer::write(*(Quaternion*)instance);
        }
        // base class
        static std::vector<ReflectionInstance> getQuaternionBaseClassReflectionInstanceList(void* instance){
            std::vector<ReflectionInstance> out_list;
            
            return out_list;
        }
        // fields
        static const char* getFieldName_w(){ return "w";}
        static const char* getFieldTypeName_w(){ return "float";}
        static void set_w(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<Quaternion*>(instance);
            typed_instance->w = *static_cast<float*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "w");
        }
        static void* get_w(void* instance){ return static_cast<void*>(&(static_cast<Quaternion*>(instance)->w));}
        static bool isArray_w(){ return false; }
        static const char* getFieldName_x(){ return "x";}
        static const char* getFieldTypeName_x(){ return "float";}
        static void set_x(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<Quaternion*>(instance);
            typed_instance->x = *static_cast<float*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "x");
        }
        static void* get_x(void* instance){ return static_cast<void*>(&(static_cast<Quaternion*>(instance)->x));}
        static bool isArray_x(){ return false; }
        static const char* getFieldName_y(){ return "y";}
        static const char* getFieldTypeName_y(){ return "float";}
        static void set_y(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<Quaternion*>(instance);
            typed_instance->y = *static_cast<float*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "y");
        }
        static void* get_y(void* instance){ return static_cast<void*>(&(static_cast<Quaternion*>(instance)->y));}
        static bool isArray_y(){ return false; }
        static const char* getFieldName_z(){ return "z";}
        static const char* getFieldTypeName_z(){ return "float";}
        static void set_z(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<Quaternion*>(instance);
            typed_instance->z = *static_cast<float*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "z");
        }
        static void* get_z(void* instance){ return static_cast<void*>(&(static_cast<Quaternion*>(instance)->z));}
        static bool isArray_z(){ return false; }
    };
}//namespace TypeFieldReflectionOparator


    void TypeWrapperRegister_Quaternion(){
        std::unordered_map<std::string, std::string> meta_tags_w;
        

        FieldFunctionTuple* f_field_function_tuple_w=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeQuaternionOperator::set_w,
            &TypeFieldReflectionOparator::TypeQuaternionOperator::get_w,
            &TypeFieldReflectionOparator::TypeQuaternionOperator::getClassName,
            &TypeFieldReflectionOparator::TypeQuaternionOperator::getFieldName_w,
            &TypeFieldReflectionOparator::TypeQuaternionOperator::getFieldTypeName_w,
            &TypeFieldReflectionOparator::TypeQuaternionOperator::isArray_w,
            meta_tags_w
        );
        REGISTER_FIELD_TO_MAP("Quaternion", f_field_function_tuple_w);
        std::unordered_map<std::string, std::string> meta_tags_x;
        

        FieldFunctionTuple* f_field_function_tuple_x=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeQuaternionOperator::set_x,
            &TypeFieldReflectionOparator::TypeQuaternionOperator::get_x,
            &TypeFieldReflectionOparator::TypeQuaternionOperator::getClassName,
            &TypeFieldReflectionOparator::TypeQuaternionOperator::getFieldName_x,
            &TypeFieldReflectionOparator::TypeQuaternionOperator::getFieldTypeName_x,
            &TypeFieldReflectionOparator::TypeQuaternionOperator::isArray_x,
            meta_tags_x
        );
        REGISTER_FIELD_TO_MAP("Quaternion", f_field_function_tuple_x);
        std::unordered_map<std::string, std::string> meta_tags_y;
        

        FieldFunctionTuple* f_field_function_tuple_y=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeQuaternionOperator::set_y,
            &TypeFieldReflectionOparator::TypeQuaternionOperator::get_y,
            &TypeFieldReflectionOparator::TypeQuaternionOperator::getClassName,
            &TypeFieldReflectionOparator::TypeQuaternionOperator::getFieldName_y,
            &TypeFieldReflectionOparator::TypeQuaternionOperator::getFieldTypeName_y,
            &TypeFieldReflectionOparator::TypeQuaternionOperator::isArray_y,
            meta_tags_y
        );
        REGISTER_FIELD_TO_MAP("Quaternion", f_field_function_tuple_y);
        std::unordered_map<std::string, std::string> meta_tags_z;
        

        FieldFunctionTuple* f_field_function_tuple_z=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeQuaternionOperator::set_z,
            &TypeFieldReflectionOparator::TypeQuaternionOperator::get_z,
            &TypeFieldReflectionOparator::TypeQuaternionOperator::getClassName,
            &TypeFieldReflectionOparator::TypeQuaternionOperator::getFieldName_z,
            &TypeFieldReflectionOparator::TypeQuaternionOperator::getFieldTypeName_z,
            &TypeFieldReflectionOparator::TypeQuaternionOperator::isArray_z,
            meta_tags_z
        );
        REGISTER_FIELD_TO_MAP("Quaternion", f_field_function_tuple_z);
        
        
        ClassFunctionTuple* f_class_function_tuple_Quaternion=new ClassFunctionTuple(
            &TypeFieldReflectionOparator::TypeQuaternionOperator::getQuaternionBaseClassReflectionInstanceList,
            &TypeFieldReflectionOparator::TypeQuaternionOperator::constructorWithJson,
            &TypeFieldReflectionOparator::TypeQuaternionOperator::writeByName);
        REGISTER_BASE_CLASS_TO_MAP("Quaternion", f_class_function_tuple_Quaternion);
    }
namespace TypeWrappersRegister{
        void Quaternion(){ TypeWrapperRegister_Quaternion();}
}//namespace TypeWrappersRegister

}//namespace Reflection
}//namespace Lizeral


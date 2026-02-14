#pragma once
#include "runtime\core\math\axis_aligned.h"

namespace Lizeral{
    class AxisAlignedBox;
namespace Reflection{
namespace TypeFieldReflectionOparator{
    class TypeAxisAlignedBoxOperator{
    public:
        static const char* getClassName(){ return "AxisAlignedBox";}
        static void* constructorWithJson(const PJson& json_context){
            AxisAlignedBox* ret_instance= new AxisAlignedBox;
            PSerializer::read(json_context, *ret_instance);
            return ret_instance;
        }
        static PJson writeByName(void* instance){
            return PSerializer::write(*(AxisAlignedBox*)instance);
        }
        // base class
        static int getAxisAlignedBoxBaseClassReflectionInstanceList(ReflectionInstance* &out_list, void* instance){
            int count = 0;
            
            return count;
        }
        // fields
        static const char* getFieldName_m_center(){ return "m_center";}
        static const char* getFieldTypeName_m_center(){ return "Lizeral::Vector3";}
        static void set_m_center(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<AxisAlignedBox*>(instance);
            typed_instance->m_center = *static_cast<Lizeral::Vector3*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "m_center");
        }
        static void* get_m_center(void* instance){ return static_cast<void*>(&(static_cast<AxisAlignedBox*>(instance)->m_center));}
        static bool isArray_m_center(){ return false; }
        static const char* getFieldName_m_half_extent(){ return "m_half_extent";}
        static const char* getFieldTypeName_m_half_extent(){ return "Lizeral::Vector3";}
        static void set_m_half_extent(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<AxisAlignedBox*>(instance);
            typed_instance->m_half_extent = *static_cast<Lizeral::Vector3*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "m_half_extent");
        }
        static void* get_m_half_extent(void* instance){ return static_cast<void*>(&(static_cast<AxisAlignedBox*>(instance)->m_half_extent));}
        static bool isArray_m_half_extent(){ return false; }
        static const char* getFieldName_m_min_corner(){ return "m_min_corner";}
        static const char* getFieldTypeName_m_min_corner(){ return "Lizeral::Vector3";}
        static void set_m_min_corner(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<AxisAlignedBox*>(instance);
            typed_instance->m_min_corner = *static_cast<Lizeral::Vector3*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "m_min_corner");
        }
        static void* get_m_min_corner(void* instance){ return static_cast<void*>(&(static_cast<AxisAlignedBox*>(instance)->m_min_corner));}
        static bool isArray_m_min_corner(){ return false; }
        static const char* getFieldName_m_max_corner(){ return "m_max_corner";}
        static const char* getFieldTypeName_m_max_corner(){ return "Lizeral::Vector3";}
        static void set_m_max_corner(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<AxisAlignedBox*>(instance);
            typed_instance->m_max_corner = *static_cast<Lizeral::Vector3*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "m_max_corner");
        }
        static void* get_m_max_corner(void* instance){ return static_cast<void*>(&(static_cast<AxisAlignedBox*>(instance)->m_max_corner));}
        static bool isArray_m_max_corner(){ return false; }
    };
}//namespace TypeFieldReflectionOparator


    void TypeWrapperRegister_AxisAlignedBox(){
        FieldFunctionTuple* f_field_function_tuple_m_center=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeAxisAlignedBoxOperator::set_m_center,
            &TypeFieldReflectionOparator::TypeAxisAlignedBoxOperator::get_m_center,
            &TypeFieldReflectionOparator::TypeAxisAlignedBoxOperator::getClassName,
            &TypeFieldReflectionOparator::TypeAxisAlignedBoxOperator::getFieldName_m_center,
            &TypeFieldReflectionOparator::TypeAxisAlignedBoxOperator::getFieldTypeName_m_center,
            &TypeFieldReflectionOparator::TypeAxisAlignedBoxOperator::isArray_m_center);
        REGISTER_FIELD_TO_MAP("AxisAlignedBox", f_field_function_tuple_m_center);
        FieldFunctionTuple* f_field_function_tuple_m_half_extent=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeAxisAlignedBoxOperator::set_m_half_extent,
            &TypeFieldReflectionOparator::TypeAxisAlignedBoxOperator::get_m_half_extent,
            &TypeFieldReflectionOparator::TypeAxisAlignedBoxOperator::getClassName,
            &TypeFieldReflectionOparator::TypeAxisAlignedBoxOperator::getFieldName_m_half_extent,
            &TypeFieldReflectionOparator::TypeAxisAlignedBoxOperator::getFieldTypeName_m_half_extent,
            &TypeFieldReflectionOparator::TypeAxisAlignedBoxOperator::isArray_m_half_extent);
        REGISTER_FIELD_TO_MAP("AxisAlignedBox", f_field_function_tuple_m_half_extent);
        FieldFunctionTuple* f_field_function_tuple_m_min_corner=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeAxisAlignedBoxOperator::set_m_min_corner,
            &TypeFieldReflectionOparator::TypeAxisAlignedBoxOperator::get_m_min_corner,
            &TypeFieldReflectionOparator::TypeAxisAlignedBoxOperator::getClassName,
            &TypeFieldReflectionOparator::TypeAxisAlignedBoxOperator::getFieldName_m_min_corner,
            &TypeFieldReflectionOparator::TypeAxisAlignedBoxOperator::getFieldTypeName_m_min_corner,
            &TypeFieldReflectionOparator::TypeAxisAlignedBoxOperator::isArray_m_min_corner);
        REGISTER_FIELD_TO_MAP("AxisAlignedBox", f_field_function_tuple_m_min_corner);
        FieldFunctionTuple* f_field_function_tuple_m_max_corner=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeAxisAlignedBoxOperator::set_m_max_corner,
            &TypeFieldReflectionOparator::TypeAxisAlignedBoxOperator::get_m_max_corner,
            &TypeFieldReflectionOparator::TypeAxisAlignedBoxOperator::getClassName,
            &TypeFieldReflectionOparator::TypeAxisAlignedBoxOperator::getFieldName_m_max_corner,
            &TypeFieldReflectionOparator::TypeAxisAlignedBoxOperator::getFieldTypeName_m_max_corner,
            &TypeFieldReflectionOparator::TypeAxisAlignedBoxOperator::isArray_m_max_corner);
        REGISTER_FIELD_TO_MAP("AxisAlignedBox", f_field_function_tuple_m_max_corner);
        
        
        ClassFunctionTuple* f_class_function_tuple_AxisAlignedBox=new ClassFunctionTuple(
            &TypeFieldReflectionOparator::TypeAxisAlignedBoxOperator::getAxisAlignedBoxBaseClassReflectionInstanceList,
            &TypeFieldReflectionOparator::TypeAxisAlignedBoxOperator::constructorWithJson,
            &TypeFieldReflectionOparator::TypeAxisAlignedBoxOperator::writeByName);
        REGISTER_BASE_CLASS_TO_MAP("AxisAlignedBox", f_class_function_tuple_AxisAlignedBox);
    }
namespace TypeWrappersRegister{
        void AxisAlignedBox(){ TypeWrapperRegister_AxisAlignedBox();}
}//namespace TypeWrappersRegister

}//namespace Reflection
}//namespace Lizeral


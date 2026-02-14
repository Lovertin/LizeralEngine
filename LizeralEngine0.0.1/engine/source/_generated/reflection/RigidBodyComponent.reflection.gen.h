#pragma once
#include "runtime\function\ecs\components\RigidBody\RigidBodyComponent.h"

namespace Lizeral{
    class RigidBodyComponent;
namespace Reflection{
namespace TypeFieldReflectionOparator{
    class TypeRigidBodyComponentOperator{
    public:
        static const char* getClassName(){ return "RigidBodyComponent";}
        static void* constructorWithJson(const PJson& json_context){
            RigidBodyComponent* ret_instance= new RigidBodyComponent;
            PSerializer::read(json_context, *ret_instance);
            return ret_instance;
        }
        static PJson writeByName(void* instance){
            return PSerializer::write(*(RigidBodyComponent*)instance);
        }
        // base class
        static int getRigidBodyComponentBaseClassReflectionInstanceList(ReflectionInstance* &out_list, void* instance){
            int count = 1;
            out_list = new ReflectionInstance[count];
            for (int i=0;i<count;++i){
               out_list[i] = TypeMetaDef(Lizeral::Component,static_cast<RigidBodyComponent*>(instance));
            }
            return count;
        }
        // fields
        static const char* getFieldName_m_mass(){ return "m_mass";}
        static const char* getFieldTypeName_m_mass(){ return "float";}
        static void set_m_mass(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<RigidBodyComponent*>(instance);
            typed_instance->m_mass = *static_cast<float*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "m_mass");
        }
        static void* get_m_mass(void* instance){ return static_cast<void*>(&(static_cast<RigidBodyComponent*>(instance)->m_mass));}
        static bool isArray_m_mass(){ return false; }
        static const char* getFieldName_m_friction(){ return "m_friction";}
        static const char* getFieldTypeName_m_friction(){ return "float";}
        static void set_m_friction(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<RigidBodyComponent*>(instance);
            typed_instance->m_friction = *static_cast<float*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "m_friction");
        }
        static void* get_m_friction(void* instance){ return static_cast<void*>(&(static_cast<RigidBodyComponent*>(instance)->m_friction));}
        static bool isArray_m_friction(){ return false; }
        static const char* getFieldName_m_is_kinematic(){ return "m_is_kinematic";}
        static const char* getFieldTypeName_m_is_kinematic(){ return "bool";}
        static void set_m_is_kinematic(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<RigidBodyComponent*>(instance);
            typed_instance->m_is_kinematic = *static_cast<bool*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "m_is_kinematic");
        }
        static void* get_m_is_kinematic(void* instance){ return static_cast<void*>(&(static_cast<RigidBodyComponent*>(instance)->m_is_kinematic));}
        static bool isArray_m_is_kinematic(){ return false; }
    };
}//namespace TypeFieldReflectionOparator


    void TypeWrapperRegister_RigidBodyComponent(){
        FieldFunctionTuple* f_field_function_tuple_m_mass=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeRigidBodyComponentOperator::set_m_mass,
            &TypeFieldReflectionOparator::TypeRigidBodyComponentOperator::get_m_mass,
            &TypeFieldReflectionOparator::TypeRigidBodyComponentOperator::getClassName,
            &TypeFieldReflectionOparator::TypeRigidBodyComponentOperator::getFieldName_m_mass,
            &TypeFieldReflectionOparator::TypeRigidBodyComponentOperator::getFieldTypeName_m_mass,
            &TypeFieldReflectionOparator::TypeRigidBodyComponentOperator::isArray_m_mass);
        REGISTER_FIELD_TO_MAP("RigidBodyComponent", f_field_function_tuple_m_mass);
        FieldFunctionTuple* f_field_function_tuple_m_friction=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeRigidBodyComponentOperator::set_m_friction,
            &TypeFieldReflectionOparator::TypeRigidBodyComponentOperator::get_m_friction,
            &TypeFieldReflectionOparator::TypeRigidBodyComponentOperator::getClassName,
            &TypeFieldReflectionOparator::TypeRigidBodyComponentOperator::getFieldName_m_friction,
            &TypeFieldReflectionOparator::TypeRigidBodyComponentOperator::getFieldTypeName_m_friction,
            &TypeFieldReflectionOparator::TypeRigidBodyComponentOperator::isArray_m_friction);
        REGISTER_FIELD_TO_MAP("RigidBodyComponent", f_field_function_tuple_m_friction);
        FieldFunctionTuple* f_field_function_tuple_m_is_kinematic=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeRigidBodyComponentOperator::set_m_is_kinematic,
            &TypeFieldReflectionOparator::TypeRigidBodyComponentOperator::get_m_is_kinematic,
            &TypeFieldReflectionOparator::TypeRigidBodyComponentOperator::getClassName,
            &TypeFieldReflectionOparator::TypeRigidBodyComponentOperator::getFieldName_m_is_kinematic,
            &TypeFieldReflectionOparator::TypeRigidBodyComponentOperator::getFieldTypeName_m_is_kinematic,
            &TypeFieldReflectionOparator::TypeRigidBodyComponentOperator::isArray_m_is_kinematic);
        REGISTER_FIELD_TO_MAP("RigidBodyComponent", f_field_function_tuple_m_is_kinematic);
        
        
        ClassFunctionTuple* f_class_function_tuple_RigidBodyComponent=new ClassFunctionTuple(
            &TypeFieldReflectionOparator::TypeRigidBodyComponentOperator::getRigidBodyComponentBaseClassReflectionInstanceList,
            &TypeFieldReflectionOparator::TypeRigidBodyComponentOperator::constructorWithJson,
            &TypeFieldReflectionOparator::TypeRigidBodyComponentOperator::writeByName);
        REGISTER_BASE_CLASS_TO_MAP("RigidBodyComponent", f_class_function_tuple_RigidBodyComponent);
    }
namespace TypeWrappersRegister{
        void RigidBodyComponent(){ TypeWrapperRegister_RigidBodyComponent();}
}//namespace TypeWrappersRegister

}//namespace Reflection
}//namespace Lizeral


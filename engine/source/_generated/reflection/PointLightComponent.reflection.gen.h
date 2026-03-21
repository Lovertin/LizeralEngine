#pragma once
#include "runtime/core/meta/serializer/serializer.h"
#include "runtime\function\ecs\components\Light\PointLightComponent.h"

namespace Lizeral{
    class PointLightComponent;
namespace Reflection{
namespace TypeFieldReflectionOparator{
    class TypePointLightComponentOperator{
    public:
        static const char* getClassName(){ return "PointLightComponent";}
        static void* constructorWithJson(const PJson& json_context){
            PointLightComponent* ret_instance= new PointLightComponent;
            PSerializer::read(json_context, *ret_instance);
            return ret_instance;
        }
        static PJson writeByName(void* instance){
            return PSerializer::write(*(PointLightComponent*)instance);
        }
        // base class
        static std::vector<ReflectionInstance> getPointLightComponentBaseClassReflectionInstanceList(void* instance){
            std::vector<ReflectionInstance> out_list;
            out_list.push_back(TypeMetaDef(Lizeral::Component, static_cast<PointLightComponent*>(instance)));
            return out_list;
        }
        // fields
        static const char* getFieldName_m_spotintensity(){ return "m_spotintensity";}
        static const char* getFieldTypeName_m_spotintensity(){ return "float";}
        static void set_m_spotintensity(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<PointLightComponent*>(instance);
            typed_instance->m_spotintensity = *static_cast<float*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "m_spotintensity");
        }
        static void* get_m_spotintensity(void* instance){ return static_cast<void*>(&(static_cast<PointLightComponent*>(instance)->m_spotintensity));}
        static bool isArray_m_spotintensity(){ return false; }
        static const char* getFieldName_m_spotlightcolor(){ return "m_spotlightcolor";}
        static const char* getFieldTypeName_m_spotlightcolor(){ return "Lizeral::Vector3";}
        static void set_m_spotlightcolor(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<PointLightComponent*>(instance);
            typed_instance->m_spotlightcolor = *static_cast<Lizeral::Vector3*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "m_spotlightcolor");
        }
        static void* get_m_spotlightcolor(void* instance){ return static_cast<void*>(&(static_cast<PointLightComponent*>(instance)->m_spotlightcolor));}
        static bool isArray_m_spotlightcolor(){ return false; }
        static const char* getFieldName_m_spotradius(){ return "m_spotradius";}
        static const char* getFieldTypeName_m_spotradius(){ return "float";}
        static void set_m_spotradius(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<PointLightComponent*>(instance);
            typed_instance->m_spotradius = *static_cast<float*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "m_spotradius");
        }
        static void* get_m_spotradius(void* instance){ return static_cast<void*>(&(static_cast<PointLightComponent*>(instance)->m_spotradius));}
        static bool isArray_m_spotradius(){ return false; }
    };
}//namespace TypeFieldReflectionOparator


    void TypeWrapperRegister_PointLightComponent(){
        std::unordered_map<std::string, std::string> meta_tags_m_spotintensity;
        meta_tags_m_spotintensity.insert({"Enable", ""});

        FieldFunctionTuple* f_field_function_tuple_m_spotintensity=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypePointLightComponentOperator::set_m_spotintensity,
            &TypeFieldReflectionOparator::TypePointLightComponentOperator::get_m_spotintensity,
            &TypeFieldReflectionOparator::TypePointLightComponentOperator::getClassName,
            &TypeFieldReflectionOparator::TypePointLightComponentOperator::getFieldName_m_spotintensity,
            &TypeFieldReflectionOparator::TypePointLightComponentOperator::getFieldTypeName_m_spotintensity,
            &TypeFieldReflectionOparator::TypePointLightComponentOperator::isArray_m_spotintensity,
            meta_tags_m_spotintensity
        );
        REGISTER_FIELD_TO_MAP("PointLightComponent", f_field_function_tuple_m_spotintensity);
        std::unordered_map<std::string, std::string> meta_tags_m_spotlightcolor;
        meta_tags_m_spotlightcolor.insert({"Enable", ""});

        FieldFunctionTuple* f_field_function_tuple_m_spotlightcolor=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypePointLightComponentOperator::set_m_spotlightcolor,
            &TypeFieldReflectionOparator::TypePointLightComponentOperator::get_m_spotlightcolor,
            &TypeFieldReflectionOparator::TypePointLightComponentOperator::getClassName,
            &TypeFieldReflectionOparator::TypePointLightComponentOperator::getFieldName_m_spotlightcolor,
            &TypeFieldReflectionOparator::TypePointLightComponentOperator::getFieldTypeName_m_spotlightcolor,
            &TypeFieldReflectionOparator::TypePointLightComponentOperator::isArray_m_spotlightcolor,
            meta_tags_m_spotlightcolor
        );
        REGISTER_FIELD_TO_MAP("PointLightComponent", f_field_function_tuple_m_spotlightcolor);
        std::unordered_map<std::string, std::string> meta_tags_m_spotradius;
        meta_tags_m_spotradius.insert({"Enable", ""});

        FieldFunctionTuple* f_field_function_tuple_m_spotradius=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypePointLightComponentOperator::set_m_spotradius,
            &TypeFieldReflectionOparator::TypePointLightComponentOperator::get_m_spotradius,
            &TypeFieldReflectionOparator::TypePointLightComponentOperator::getClassName,
            &TypeFieldReflectionOparator::TypePointLightComponentOperator::getFieldName_m_spotradius,
            &TypeFieldReflectionOparator::TypePointLightComponentOperator::getFieldTypeName_m_spotradius,
            &TypeFieldReflectionOparator::TypePointLightComponentOperator::isArray_m_spotradius,
            meta_tags_m_spotradius
        );
        REGISTER_FIELD_TO_MAP("PointLightComponent", f_field_function_tuple_m_spotradius);
        
        
        ClassFunctionTuple* f_class_function_tuple_PointLightComponent=new ClassFunctionTuple(
            &TypeFieldReflectionOparator::TypePointLightComponentOperator::getPointLightComponentBaseClassReflectionInstanceList,
            &TypeFieldReflectionOparator::TypePointLightComponentOperator::constructorWithJson,
            &TypeFieldReflectionOparator::TypePointLightComponentOperator::writeByName);
        REGISTER_BASE_CLASS_TO_MAP("PointLightComponent", f_class_function_tuple_PointLightComponent);
    }
namespace TypeWrappersRegister{
        void PointLightComponent(){ TypeWrapperRegister_PointLightComponent();}
}//namespace TypeWrappersRegister

}//namespace Reflection
}//namespace Lizeral


#pragma once
#include "runtime/core/meta/serializer/serializer.h"
#include "runtime\function\ecs\components\Model\VulkanModelComponent.h"

namespace Lizeral{
    class VulkanModelComponent;
namespace Reflection{
namespace TypeFieldReflectionOparator{
    class TypeVulkanModelComponentOperator{
    public:
        static const char* getClassName(){ return "VulkanModelComponent";}
        static void* constructorWithJson(const PJson& json_context){
            VulkanModelComponent* ret_instance= new VulkanModelComponent;
            PSerializer::read(json_context, *ret_instance);
            return ret_instance;
        }
        static PJson writeByName(void* instance){
            return PSerializer::write(*(VulkanModelComponent*)instance);
        }
        // base class
        static std::vector<ReflectionInstance> getVulkanModelComponentBaseClassReflectionInstanceList(void* instance){
            std::vector<ReflectionInstance> out_list;
            out_list.push_back(TypeMetaDef(Lizeral::Component, static_cast<VulkanModelComponent*>(instance)));
            return out_list;
        }
        // fields
        static const char* getFieldName_m_ModelPath(){ return "m_ModelPath";}
        static const char* getFieldTypeName_m_ModelPath(){ return "std::string";}
        static void set_m_ModelPath(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<VulkanModelComponent*>(instance);
            typed_instance->m_ModelPath = *static_cast<std::string*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "m_ModelPath");
        }
        static void* get_m_ModelPath(void* instance){ return static_cast<void*>(&(static_cast<VulkanModelComponent*>(instance)->m_ModelPath));}
        static bool isArray_m_ModelPath(){ return false; }
    };
}//namespace TypeFieldReflectionOparator


    void TypeWrapperRegister_VulkanModelComponent(){
        std::unordered_map<std::string, std::string> meta_tags_m_ModelPath;
        meta_tags_m_ModelPath.insert({"UI", "Address"});
        meta_tags_m_ModelPath.insert({"Enable", ""});

        FieldFunctionTuple* f_field_function_tuple_m_ModelPath=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeVulkanModelComponentOperator::set_m_ModelPath,
            &TypeFieldReflectionOparator::TypeVulkanModelComponentOperator::get_m_ModelPath,
            &TypeFieldReflectionOparator::TypeVulkanModelComponentOperator::getClassName,
            &TypeFieldReflectionOparator::TypeVulkanModelComponentOperator::getFieldName_m_ModelPath,
            &TypeFieldReflectionOparator::TypeVulkanModelComponentOperator::getFieldTypeName_m_ModelPath,
            &TypeFieldReflectionOparator::TypeVulkanModelComponentOperator::isArray_m_ModelPath,
            meta_tags_m_ModelPath
        );
        REGISTER_FIELD_TO_MAP("VulkanModelComponent", f_field_function_tuple_m_ModelPath);
        
        
        ClassFunctionTuple* f_class_function_tuple_VulkanModelComponent=new ClassFunctionTuple(
            &TypeFieldReflectionOparator::TypeVulkanModelComponentOperator::getVulkanModelComponentBaseClassReflectionInstanceList,
            &TypeFieldReflectionOparator::TypeVulkanModelComponentOperator::constructorWithJson,
            &TypeFieldReflectionOparator::TypeVulkanModelComponentOperator::writeByName);
        REGISTER_BASE_CLASS_TO_MAP("VulkanModelComponent", f_class_function_tuple_VulkanModelComponent);
    }
namespace TypeWrappersRegister{
        void VulkanModelComponent(){ TypeWrapperRegister_VulkanModelComponent();}
}//namespace TypeWrappersRegister

}//namespace Reflection
}//namespace Lizeral


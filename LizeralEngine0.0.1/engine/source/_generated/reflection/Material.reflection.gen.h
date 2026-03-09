#pragma once
#include "runtime\function\res_type\Material\Material.h"

namespace Lizeral{
    class Material;
namespace Reflection{
namespace TypeFieldReflectionOparator{
    class TypeMaterialOperator{
    public:
        static const char* getClassName(){ return "Material";}
        static void* constructorWithJson(const PJson& json_context){
            Material* ret_instance= new Material;
            PSerializer::read(json_context, *ret_instance);
            return ret_instance;
        }
        static PJson writeByName(void* instance){
            return PSerializer::write(*(Material*)instance);
        }
        // base class
        static std::vector<ReflectionInstance> getMaterialBaseClassReflectionInstanceList(void* instance){
            std::vector<ReflectionInstance> out_list;
            
            return out_list;
        }
        // fields
        static const char* getFieldName_m_VertShaderPath(){ return "m_VertShaderPath";}
        static const char* getFieldTypeName_m_VertShaderPath(){ return "std::string";}
        static void set_m_VertShaderPath(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<Material*>(instance);
            typed_instance->m_VertShaderPath = *static_cast<std::string*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "m_VertShaderPath");
        }
        static void* get_m_VertShaderPath(void* instance){ return static_cast<void*>(&(static_cast<Material*>(instance)->m_VertShaderPath));}
        static bool isArray_m_VertShaderPath(){ return false; }
        static const char* getFieldName_m_FragShaderPath(){ return "m_FragShaderPath";}
        static const char* getFieldTypeName_m_FragShaderPath(){ return "std::string";}
        static void set_m_FragShaderPath(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<Material*>(instance);
            typed_instance->m_FragShaderPath = *static_cast<std::string*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "m_FragShaderPath");
        }
        static void* get_m_FragShaderPath(void* instance){ return static_cast<void*>(&(static_cast<Material*>(instance)->m_FragShaderPath));}
        static bool isArray_m_FragShaderPath(){ return false; }
    };
}//namespace TypeFieldReflectionOparator


    void TypeWrapperRegister_Material(){
        std::unordered_map<std::string, std::string> meta_tags_m_VertShaderPath;
        meta_tags_m_VertShaderPath.insert({"Enable", ""});

        FieldFunctionTuple* f_field_function_tuple_m_VertShaderPath=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeMaterialOperator::set_m_VertShaderPath,
            &TypeFieldReflectionOparator::TypeMaterialOperator::get_m_VertShaderPath,
            &TypeFieldReflectionOparator::TypeMaterialOperator::getClassName,
            &TypeFieldReflectionOparator::TypeMaterialOperator::getFieldName_m_VertShaderPath,
            &TypeFieldReflectionOparator::TypeMaterialOperator::getFieldTypeName_m_VertShaderPath,
            &TypeFieldReflectionOparator::TypeMaterialOperator::isArray_m_VertShaderPath,
            meta_tags_m_VertShaderPath
        );
        REGISTER_FIELD_TO_MAP("Material", f_field_function_tuple_m_VertShaderPath);
        std::unordered_map<std::string, std::string> meta_tags_m_FragShaderPath;
        meta_tags_m_FragShaderPath.insert({"Enable", ""});

        FieldFunctionTuple* f_field_function_tuple_m_FragShaderPath=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeMaterialOperator::set_m_FragShaderPath,
            &TypeFieldReflectionOparator::TypeMaterialOperator::get_m_FragShaderPath,
            &TypeFieldReflectionOparator::TypeMaterialOperator::getClassName,
            &TypeFieldReflectionOparator::TypeMaterialOperator::getFieldName_m_FragShaderPath,
            &TypeFieldReflectionOparator::TypeMaterialOperator::getFieldTypeName_m_FragShaderPath,
            &TypeFieldReflectionOparator::TypeMaterialOperator::isArray_m_FragShaderPath,
            meta_tags_m_FragShaderPath
        );
        REGISTER_FIELD_TO_MAP("Material", f_field_function_tuple_m_FragShaderPath);
        
        
        ClassFunctionTuple* f_class_function_tuple_Material=new ClassFunctionTuple(
            &TypeFieldReflectionOparator::TypeMaterialOperator::getMaterialBaseClassReflectionInstanceList,
            &TypeFieldReflectionOparator::TypeMaterialOperator::constructorWithJson,
            &TypeFieldReflectionOparator::TypeMaterialOperator::writeByName);
        REGISTER_BASE_CLASS_TO_MAP("Material", f_class_function_tuple_Material);
    }
namespace TypeWrappersRegister{
        void Material(){ TypeWrapperRegister_Material();}
}//namespace TypeWrappersRegister

}//namespace Reflection
}//namespace Lizeral


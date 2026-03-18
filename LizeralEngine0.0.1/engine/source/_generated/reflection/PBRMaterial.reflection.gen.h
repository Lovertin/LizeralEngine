#pragma once
#include "runtime/core/meta/serializer/serializer.h"
#include "runtime\function\res_type\Material\PBRMaterial.h"

namespace Lizeral{
    class PBRMaterial;
namespace Reflection{
namespace TypeFieldReflectionOparator{
    class TypePBRMaterialOperator{
    public:
        static const char* getClassName(){ return "PBRMaterial";}
        static void* constructorWithJson(const PJson& json_context){
            PBRMaterial* ret_instance= new PBRMaterial;
            PSerializer::read(json_context, *ret_instance);
            return ret_instance;
        }
        static PJson writeByName(void* instance){
            return PSerializer::write(*(PBRMaterial*)instance);
        }
        // base class
        static std::vector<ReflectionInstance> getPBRMaterialBaseClassReflectionInstanceList(void* instance){
            std::vector<ReflectionInstance> out_list;
            out_list.push_back(TypeMetaDef(Lizeral::Material, static_cast<PBRMaterial*>(instance)));
            return out_list;
        }
        // fields
        static const char* getFieldName_m_Albedo(){ return "m_Albedo";}
        static const char* getFieldTypeName_m_Albedo(){ return "Lizeral::Vector3";}
        static void set_m_Albedo(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<PBRMaterial*>(instance);
            typed_instance->m_Albedo = *static_cast<Lizeral::Vector3*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "m_Albedo");
        }
        static void* get_m_Albedo(void* instance){ return static_cast<void*>(&(static_cast<PBRMaterial*>(instance)->m_Albedo));}
        static bool isArray_m_Albedo(){ return false; }
        static const char* getFieldName_m_Metallic(){ return "m_Metallic";}
        static const char* getFieldTypeName_m_Metallic(){ return "float";}
        static void set_m_Metallic(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<PBRMaterial*>(instance);
            typed_instance->m_Metallic = *static_cast<float*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "m_Metallic");
        }
        static void* get_m_Metallic(void* instance){ return static_cast<void*>(&(static_cast<PBRMaterial*>(instance)->m_Metallic));}
        static bool isArray_m_Metallic(){ return false; }
        static const char* getFieldName_m_Roughness(){ return "m_Roughness";}
        static const char* getFieldTypeName_m_Roughness(){ return "float";}
        static void set_m_Roughness(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<PBRMaterial*>(instance);
            typed_instance->m_Roughness = *static_cast<float*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "m_Roughness");
        }
        static void* get_m_Roughness(void* instance){ return static_cast<void*>(&(static_cast<PBRMaterial*>(instance)->m_Roughness));}
        static bool isArray_m_Roughness(){ return false; }
        static const char* getFieldName_m_AO(){ return "m_AO";}
        static const char* getFieldTypeName_m_AO(){ return "float";}
        static void set_m_AO(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<PBRMaterial*>(instance);
            typed_instance->m_AO = *static_cast<float*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "m_AO");
        }
        static void* get_m_AO(void* instance){ return static_cast<void*>(&(static_cast<PBRMaterial*>(instance)->m_AO));}
        static bool isArray_m_AO(){ return false; }
        static const char* getFieldName_m_AlbedoMapPath(){ return "m_AlbedoMapPath";}
        static const char* getFieldTypeName_m_AlbedoMapPath(){ return "std::string";}
        static void set_m_AlbedoMapPath(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<PBRMaterial*>(instance);
            typed_instance->m_AlbedoMapPath = *static_cast<std::string*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "m_AlbedoMapPath");
        }
        static void* get_m_AlbedoMapPath(void* instance){ return static_cast<void*>(&(static_cast<PBRMaterial*>(instance)->m_AlbedoMapPath));}
        static bool isArray_m_AlbedoMapPath(){ return false; }
        static const char* getFieldName_m_IrradianceMapPath(){ return "m_IrradianceMapPath";}
        static const char* getFieldTypeName_m_IrradianceMapPath(){ return "std::string";}
        static void set_m_IrradianceMapPath(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<PBRMaterial*>(instance);
            typed_instance->m_IrradianceMapPath = *static_cast<std::string*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "m_IrradianceMapPath");
        }
        static void* get_m_IrradianceMapPath(void* instance){ return static_cast<void*>(&(static_cast<PBRMaterial*>(instance)->m_IrradianceMapPath));}
        static bool isArray_m_IrradianceMapPath(){ return false; }
        static const char* getFieldName_m_PrefilterMapPath(){ return "m_PrefilterMapPath";}
        static const char* getFieldTypeName_m_PrefilterMapPath(){ return "std::string";}
        static void set_m_PrefilterMapPath(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<PBRMaterial*>(instance);
            typed_instance->m_PrefilterMapPath = *static_cast<std::string*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "m_PrefilterMapPath");
        }
        static void* get_m_PrefilterMapPath(void* instance){ return static_cast<void*>(&(static_cast<PBRMaterial*>(instance)->m_PrefilterMapPath));}
        static bool isArray_m_PrefilterMapPath(){ return false; }
    };
}//namespace TypeFieldReflectionOparator


    void TypeWrapperRegister_PBRMaterial(){
        std::unordered_map<std::string, std::string> meta_tags_m_Albedo;
        meta_tags_m_Albedo.insert({"Enable", ""});

        FieldFunctionTuple* f_field_function_tuple_m_Albedo=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypePBRMaterialOperator::set_m_Albedo,
            &TypeFieldReflectionOparator::TypePBRMaterialOperator::get_m_Albedo,
            &TypeFieldReflectionOparator::TypePBRMaterialOperator::getClassName,
            &TypeFieldReflectionOparator::TypePBRMaterialOperator::getFieldName_m_Albedo,
            &TypeFieldReflectionOparator::TypePBRMaterialOperator::getFieldTypeName_m_Albedo,
            &TypeFieldReflectionOparator::TypePBRMaterialOperator::isArray_m_Albedo,
            meta_tags_m_Albedo
        );
        REGISTER_FIELD_TO_MAP("PBRMaterial", f_field_function_tuple_m_Albedo);
        std::unordered_map<std::string, std::string> meta_tags_m_Metallic;
        meta_tags_m_Metallic.insert({"Enable", ""});

        FieldFunctionTuple* f_field_function_tuple_m_Metallic=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypePBRMaterialOperator::set_m_Metallic,
            &TypeFieldReflectionOparator::TypePBRMaterialOperator::get_m_Metallic,
            &TypeFieldReflectionOparator::TypePBRMaterialOperator::getClassName,
            &TypeFieldReflectionOparator::TypePBRMaterialOperator::getFieldName_m_Metallic,
            &TypeFieldReflectionOparator::TypePBRMaterialOperator::getFieldTypeName_m_Metallic,
            &TypeFieldReflectionOparator::TypePBRMaterialOperator::isArray_m_Metallic,
            meta_tags_m_Metallic
        );
        REGISTER_FIELD_TO_MAP("PBRMaterial", f_field_function_tuple_m_Metallic);
        std::unordered_map<std::string, std::string> meta_tags_m_Roughness;
        meta_tags_m_Roughness.insert({"Enable", ""});

        FieldFunctionTuple* f_field_function_tuple_m_Roughness=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypePBRMaterialOperator::set_m_Roughness,
            &TypeFieldReflectionOparator::TypePBRMaterialOperator::get_m_Roughness,
            &TypeFieldReflectionOparator::TypePBRMaterialOperator::getClassName,
            &TypeFieldReflectionOparator::TypePBRMaterialOperator::getFieldName_m_Roughness,
            &TypeFieldReflectionOparator::TypePBRMaterialOperator::getFieldTypeName_m_Roughness,
            &TypeFieldReflectionOparator::TypePBRMaterialOperator::isArray_m_Roughness,
            meta_tags_m_Roughness
        );
        REGISTER_FIELD_TO_MAP("PBRMaterial", f_field_function_tuple_m_Roughness);
        std::unordered_map<std::string, std::string> meta_tags_m_AO;
        meta_tags_m_AO.insert({"Enable", ""});

        FieldFunctionTuple* f_field_function_tuple_m_AO=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypePBRMaterialOperator::set_m_AO,
            &TypeFieldReflectionOparator::TypePBRMaterialOperator::get_m_AO,
            &TypeFieldReflectionOparator::TypePBRMaterialOperator::getClassName,
            &TypeFieldReflectionOparator::TypePBRMaterialOperator::getFieldName_m_AO,
            &TypeFieldReflectionOparator::TypePBRMaterialOperator::getFieldTypeName_m_AO,
            &TypeFieldReflectionOparator::TypePBRMaterialOperator::isArray_m_AO,
            meta_tags_m_AO
        );
        REGISTER_FIELD_TO_MAP("PBRMaterial", f_field_function_tuple_m_AO);
        std::unordered_map<std::string, std::string> meta_tags_m_AlbedoMapPath;
        meta_tags_m_AlbedoMapPath.insert({"Enable", ""});

        FieldFunctionTuple* f_field_function_tuple_m_AlbedoMapPath=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypePBRMaterialOperator::set_m_AlbedoMapPath,
            &TypeFieldReflectionOparator::TypePBRMaterialOperator::get_m_AlbedoMapPath,
            &TypeFieldReflectionOparator::TypePBRMaterialOperator::getClassName,
            &TypeFieldReflectionOparator::TypePBRMaterialOperator::getFieldName_m_AlbedoMapPath,
            &TypeFieldReflectionOparator::TypePBRMaterialOperator::getFieldTypeName_m_AlbedoMapPath,
            &TypeFieldReflectionOparator::TypePBRMaterialOperator::isArray_m_AlbedoMapPath,
            meta_tags_m_AlbedoMapPath
        );
        REGISTER_FIELD_TO_MAP("PBRMaterial", f_field_function_tuple_m_AlbedoMapPath);
        std::unordered_map<std::string, std::string> meta_tags_m_IrradianceMapPath;
        meta_tags_m_IrradianceMapPath.insert({"Enable", ""});

        FieldFunctionTuple* f_field_function_tuple_m_IrradianceMapPath=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypePBRMaterialOperator::set_m_IrradianceMapPath,
            &TypeFieldReflectionOparator::TypePBRMaterialOperator::get_m_IrradianceMapPath,
            &TypeFieldReflectionOparator::TypePBRMaterialOperator::getClassName,
            &TypeFieldReflectionOparator::TypePBRMaterialOperator::getFieldName_m_IrradianceMapPath,
            &TypeFieldReflectionOparator::TypePBRMaterialOperator::getFieldTypeName_m_IrradianceMapPath,
            &TypeFieldReflectionOparator::TypePBRMaterialOperator::isArray_m_IrradianceMapPath,
            meta_tags_m_IrradianceMapPath
        );
        REGISTER_FIELD_TO_MAP("PBRMaterial", f_field_function_tuple_m_IrradianceMapPath);
        std::unordered_map<std::string, std::string> meta_tags_m_PrefilterMapPath;
        meta_tags_m_PrefilterMapPath.insert({"Enable", ""});

        FieldFunctionTuple* f_field_function_tuple_m_PrefilterMapPath=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypePBRMaterialOperator::set_m_PrefilterMapPath,
            &TypeFieldReflectionOparator::TypePBRMaterialOperator::get_m_PrefilterMapPath,
            &TypeFieldReflectionOparator::TypePBRMaterialOperator::getClassName,
            &TypeFieldReflectionOparator::TypePBRMaterialOperator::getFieldName_m_PrefilterMapPath,
            &TypeFieldReflectionOparator::TypePBRMaterialOperator::getFieldTypeName_m_PrefilterMapPath,
            &TypeFieldReflectionOparator::TypePBRMaterialOperator::isArray_m_PrefilterMapPath,
            meta_tags_m_PrefilterMapPath
        );
        REGISTER_FIELD_TO_MAP("PBRMaterial", f_field_function_tuple_m_PrefilterMapPath);
        
        
        ClassFunctionTuple* f_class_function_tuple_PBRMaterial=new ClassFunctionTuple(
            &TypeFieldReflectionOparator::TypePBRMaterialOperator::getPBRMaterialBaseClassReflectionInstanceList,
            &TypeFieldReflectionOparator::TypePBRMaterialOperator::constructorWithJson,
            &TypeFieldReflectionOparator::TypePBRMaterialOperator::writeByName);
        REGISTER_BASE_CLASS_TO_MAP("PBRMaterial", f_class_function_tuple_PBRMaterial);
    }
namespace TypeWrappersRegister{
        void PBRMaterial(){ TypeWrapperRegister_PBRMaterial();}
}//namespace TypeWrappersRegister

}//namespace Reflection
}//namespace Lizeral


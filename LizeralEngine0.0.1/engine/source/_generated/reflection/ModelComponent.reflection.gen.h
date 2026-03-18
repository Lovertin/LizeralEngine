#pragma once
#include "runtime/core/meta/serializer/serializer.h"
#include "runtime\function\ecs\components\Model\ModelComponent.h"

namespace Lizeral{
    class ModelComponent;
namespace Reflection{
namespace TypeFieldReflectionOparator{
    class TypeModelComponentOperator{
    public:
        static const char* getClassName(){ return "ModelComponent";}
        static void* constructorWithJson(const PJson& json_context){
            ModelComponent* ret_instance= new ModelComponent;
            PSerializer::read(json_context, *ret_instance);
            return ret_instance;
        }
        static PJson writeByName(void* instance){
            return PSerializer::write(*(ModelComponent*)instance);
        }
        // base class
        static std::vector<ReflectionInstance> getModelComponentBaseClassReflectionInstanceList(void* instance){
            std::vector<ReflectionInstance> out_list;
            out_list.push_back(TypeMetaDef(Lizeral::Component, static_cast<ModelComponent*>(instance)));
            return out_list;
        }
        // fields
        static const char* getFieldName_m_ModelPath(){ return "m_ModelPath";}
        static const char* getFieldTypeName_m_ModelPath(){ return "std::string";}
        static void set_m_ModelPath(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<ModelComponent*>(instance);
            typed_instance->m_ModelPath = *static_cast<std::string*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "m_ModelPath");
        }
        static void* get_m_ModelPath(void* instance){ return static_cast<void*>(&(static_cast<ModelComponent*>(instance)->m_ModelPath));}
        static bool isArray_m_ModelPath(){ return false; }
        static const char* getFieldName_m_OverrideMaterialPaths(){ return "m_OverrideMaterialPaths";}
        static const char* getFieldTypeName_m_OverrideMaterialPaths(){ return "std::vector<std::string>";}
        static void set_m_OverrideMaterialPaths(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<ModelComponent*>(instance);
            typed_instance->m_OverrideMaterialPaths = *static_cast<std::vector<std::string>*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "m_OverrideMaterialPaths");
        }
        static void* get_m_OverrideMaterialPaths(void* instance){ return static_cast<void*>(&(static_cast<ModelComponent*>(instance)->m_OverrideMaterialPaths));}
        static bool isArray_m_OverrideMaterialPaths(){ return true; }
        static const char* getFieldName_m_UseGlobalMaterial(){ return "m_UseGlobalMaterial";}
        static const char* getFieldTypeName_m_UseGlobalMaterial(){ return "bool";}
        static void set_m_UseGlobalMaterial(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<ModelComponent*>(instance);
            typed_instance->m_UseGlobalMaterial = *static_cast<bool*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "m_UseGlobalMaterial");
        }
        static void* get_m_UseGlobalMaterial(void* instance){ return static_cast<void*>(&(static_cast<ModelComponent*>(instance)->m_UseGlobalMaterial));}
        static bool isArray_m_UseGlobalMaterial(){ return false; }
    };
}//namespace TypeFieldReflectionOparator
namespace ArrayReflectionOperator{
#ifndef ArraystdSSvectorLstdSSstringROperatorMACRO
#define ArraystdSSvectorLstdSSstringROperatorMACRO
    class ArraystdSSvectorLstdSSstringROperator{
        public:
            static const char* getArrayTypeName(){ return "std::vector<std::string>";}
            static const char* getElementTypeName(){ return "std::string";}
            static int getSize(void* instance){
                //todo: should check validation
                return static_cast<int>(static_cast<std::vector<std::string>*>(instance)->size());
            }
            static void* get(int index,void* instance){
                //todo: should check validation
                return static_cast<void*>(&((*static_cast<std::vector<std::string>*>(instance))[index]));
            }
            static void set(int index, void* instance, void* element_value){
                //todo: should check validation
                (*static_cast<std::vector<std::string>*>(instance))[index] = *static_cast<std::string*>(element_value);
            }
    };
#endif //ArraystdSSvectorLstdSSstringROperator
}//namespace ArrayReflectionOperator

    void TypeWrapperRegister_ModelComponent(){
        std::unordered_map<std::string, std::string> meta_tags_m_ModelPath;
        meta_tags_m_ModelPath.insert({"UI", "Address"});
        meta_tags_m_ModelPath.insert({"Enable", ""});

        FieldFunctionTuple* f_field_function_tuple_m_ModelPath=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeModelComponentOperator::set_m_ModelPath,
            &TypeFieldReflectionOparator::TypeModelComponentOperator::get_m_ModelPath,
            &TypeFieldReflectionOparator::TypeModelComponentOperator::getClassName,
            &TypeFieldReflectionOparator::TypeModelComponentOperator::getFieldName_m_ModelPath,
            &TypeFieldReflectionOparator::TypeModelComponentOperator::getFieldTypeName_m_ModelPath,
            &TypeFieldReflectionOparator::TypeModelComponentOperator::isArray_m_ModelPath,
            meta_tags_m_ModelPath
        );
        REGISTER_FIELD_TO_MAP("ModelComponent", f_field_function_tuple_m_ModelPath);
        std::unordered_map<std::string, std::string> meta_tags_m_OverrideMaterialPaths;
        meta_tags_m_OverrideMaterialPaths.insert({"UI", "Address"});
        meta_tags_m_OverrideMaterialPaths.insert({"Enable", ""});

        FieldFunctionTuple* f_field_function_tuple_m_OverrideMaterialPaths=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeModelComponentOperator::set_m_OverrideMaterialPaths,
            &TypeFieldReflectionOparator::TypeModelComponentOperator::get_m_OverrideMaterialPaths,
            &TypeFieldReflectionOparator::TypeModelComponentOperator::getClassName,
            &TypeFieldReflectionOparator::TypeModelComponentOperator::getFieldName_m_OverrideMaterialPaths,
            &TypeFieldReflectionOparator::TypeModelComponentOperator::getFieldTypeName_m_OverrideMaterialPaths,
            &TypeFieldReflectionOparator::TypeModelComponentOperator::isArray_m_OverrideMaterialPaths,
            meta_tags_m_OverrideMaterialPaths
        );
        REGISTER_FIELD_TO_MAP("ModelComponent", f_field_function_tuple_m_OverrideMaterialPaths);
        std::unordered_map<std::string, std::string> meta_tags_m_UseGlobalMaterial;
        meta_tags_m_UseGlobalMaterial.insert({"Enable", ""});

        FieldFunctionTuple* f_field_function_tuple_m_UseGlobalMaterial=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeModelComponentOperator::set_m_UseGlobalMaterial,
            &TypeFieldReflectionOparator::TypeModelComponentOperator::get_m_UseGlobalMaterial,
            &TypeFieldReflectionOparator::TypeModelComponentOperator::getClassName,
            &TypeFieldReflectionOparator::TypeModelComponentOperator::getFieldName_m_UseGlobalMaterial,
            &TypeFieldReflectionOparator::TypeModelComponentOperator::getFieldTypeName_m_UseGlobalMaterial,
            &TypeFieldReflectionOparator::TypeModelComponentOperator::isArray_m_UseGlobalMaterial,
            meta_tags_m_UseGlobalMaterial
        );
        REGISTER_FIELD_TO_MAP("ModelComponent", f_field_function_tuple_m_UseGlobalMaterial);
        
        ArrayFunctionTuple* f_array_tuple_stdSSvectorLstdSSstringR = new  ArrayFunctionTuple(
            &ArrayReflectionOperator::ArraystdSSvectorLstdSSstringROperator::set,
            &ArrayReflectionOperator::ArraystdSSvectorLstdSSstringROperator::get,
            &ArrayReflectionOperator::ArraystdSSvectorLstdSSstringROperator::getSize,
            &ArrayReflectionOperator::ArraystdSSvectorLstdSSstringROperator::getArrayTypeName,
            &ArrayReflectionOperator::ArraystdSSvectorLstdSSstringROperator::getElementTypeName);
        REGISTER_ARRAY_TO_MAP("std::vector<std::string>", f_array_tuple_stdSSvectorLstdSSstringR);
        ClassFunctionTuple* f_class_function_tuple_ModelComponent=new ClassFunctionTuple(
            &TypeFieldReflectionOparator::TypeModelComponentOperator::getModelComponentBaseClassReflectionInstanceList,
            &TypeFieldReflectionOparator::TypeModelComponentOperator::constructorWithJson,
            &TypeFieldReflectionOparator::TypeModelComponentOperator::writeByName);
        REGISTER_BASE_CLASS_TO_MAP("ModelComponent", f_class_function_tuple_ModelComponent);
    }
namespace TypeWrappersRegister{
        void ModelComponent(){ TypeWrapperRegister_ModelComponent();}
}//namespace TypeWrappersRegister

}//namespace Reflection
}//namespace Lizeral


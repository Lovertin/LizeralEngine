#pragma once
#include "runtime/core/meta/serializer/serializer.h"
#include "runtime\function\res_type\texture\Texture2D.h"

namespace Lizeral{
    class Texture2D;
namespace Reflection{
namespace TypeFieldReflectionOparator{
    class TypeTexture2DOperator{
    public:
        static const char* getClassName(){ return "Texture2D";}
        static void* constructorWithJson(const PJson& json_context){
            Texture2D* ret_instance= new Texture2D;
            PSerializer::read(json_context, *ret_instance);
            return ret_instance;
        }
        static PJson writeByName(void* instance){
            return PSerializer::write(*(Texture2D*)instance);
        }
        // base class
        static std::vector<ReflectionInstance> getTexture2DBaseClassReflectionInstanceList(void* instance){
            std::vector<ReflectionInstance> out_list;
            out_list.push_back(TypeMetaDef(Lizeral::Resource, static_cast<Texture2D*>(instance)));
            return out_list;
        }
        // fields
        static const char* getFieldName_m_Path(){ return "m_Path";}
        static const char* getFieldTypeName_m_Path(){ return "std::string";}
        static void set_m_Path(void* instance, void* field_value){ 
            auto* typed_instance = static_cast<Texture2D*>(instance);
            typed_instance->m_Path = *static_cast<std::string*>(field_value);
            Lizeral::Reflection::TryNotifyReflectionUpdated(typed_instance, "m_Path");
        }
        static void* get_m_Path(void* instance){ return static_cast<void*>(&(static_cast<Texture2D*>(instance)->m_Path));}
        static bool isArray_m_Path(){ return false; }
    };
}//namespace TypeFieldReflectionOparator


    void TypeWrapperRegister_Texture2D(){
        std::unordered_map<std::string, std::string> meta_tags_m_Path;
        meta_tags_m_Path.insert({"Enable", ""});

        FieldFunctionTuple* f_field_function_tuple_m_Path=new FieldFunctionTuple(
            &TypeFieldReflectionOparator::TypeTexture2DOperator::set_m_Path,
            &TypeFieldReflectionOparator::TypeTexture2DOperator::get_m_Path,
            &TypeFieldReflectionOparator::TypeTexture2DOperator::getClassName,
            &TypeFieldReflectionOparator::TypeTexture2DOperator::getFieldName_m_Path,
            &TypeFieldReflectionOparator::TypeTexture2DOperator::getFieldTypeName_m_Path,
            &TypeFieldReflectionOparator::TypeTexture2DOperator::isArray_m_Path,
            meta_tags_m_Path
        );
        REGISTER_FIELD_TO_MAP("Texture2D", f_field_function_tuple_m_Path);
        
        
        ClassFunctionTuple* f_class_function_tuple_Texture2D=new ClassFunctionTuple(
            &TypeFieldReflectionOparator::TypeTexture2DOperator::getTexture2DBaseClassReflectionInstanceList,
            &TypeFieldReflectionOparator::TypeTexture2DOperator::constructorWithJson,
            &TypeFieldReflectionOparator::TypeTexture2DOperator::writeByName);
        REGISTER_BASE_CLASS_TO_MAP("Texture2D", f_class_function_tuple_Texture2D);
    }
namespace TypeWrappersRegister{
        void Texture2D(){ TypeWrapperRegister_Texture2D();}
}//namespace TypeWrappersRegister

}//namespace Reflection
}//namespace Lizeral


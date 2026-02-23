#pragma once
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
        static int getPBRMaterialBaseClassReflectionInstanceList(ReflectionInstance* &out_list, void* instance){
            int count = 1;
            out_list = new ReflectionInstance[count];
            for (int i=0;i<count;++i){
               out_list[i] = TypeMetaDef(Lizeral::Material,static_cast<PBRMaterial*>(instance));
            }
            return count;
        }
        // fields
        
    };
}//namespace TypeFieldReflectionOparator


    void TypeWrapperRegister_PBRMaterial(){
        
        
        
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


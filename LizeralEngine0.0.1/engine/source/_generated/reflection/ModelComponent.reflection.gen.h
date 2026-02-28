#pragma once
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
        static int getModelComponentBaseClassReflectionInstanceList(ReflectionInstance* &out_list, void* instance){
            int count = 1;
            out_list = new ReflectionInstance[count];
            for (int i=0;i<count;++i){
               out_list[i] = TypeMetaDef(Lizeral::Component,static_cast<ModelComponent*>(instance));
            }
            return count;
        }
        // fields
        
    };
}//namespace TypeFieldReflectionOparator


    void TypeWrapperRegister_ModelComponent(){
        
        
        
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


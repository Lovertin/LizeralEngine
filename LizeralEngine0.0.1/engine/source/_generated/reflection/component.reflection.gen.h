#pragma once
#include "runtime\function\ecs\components\component.h"

namespace Lizeral{
    class Component;
namespace Reflection{
namespace TypeFieldReflectionOparator{
    class TypeComponentOperator{
    public:
        static const char* getClassName(){ return "Component";}
        static void* constructorWithJson(const PJson& json_context){
            Component* ret_instance= new Component;
            PSerializer::read(json_context, *ret_instance);
            return ret_instance;
        }
        static PJson writeByName(void* instance){
            return PSerializer::write(*(Component*)instance);
        }
        // base class
        static int getComponentBaseClassReflectionInstanceList(ReflectionInstance* &out_list, void* instance){
            int count = 0;
            
            return count;
        }
        // fields
        
    };
}//namespace TypeFieldReflectionOparator


    void TypeWrapperRegister_Component(){
        
        
        
        ClassFunctionTuple* f_class_function_tuple_Component=new ClassFunctionTuple(
            &TypeFieldReflectionOparator::TypeComponentOperator::getComponentBaseClassReflectionInstanceList,
            &TypeFieldReflectionOparator::TypeComponentOperator::constructorWithJson,
            &TypeFieldReflectionOparator::TypeComponentOperator::writeByName);
        REGISTER_BASE_CLASS_TO_MAP("Component", f_class_function_tuple_Component);
    }
namespace TypeWrappersRegister{
        void Component(){ TypeWrapperRegister_Component();}
}//namespace TypeWrappersRegister

}//namespace Reflection
}//namespace Lizeral


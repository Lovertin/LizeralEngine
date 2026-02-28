#pragma once
#include "runtime\function\res_type\Model\Mesh.h"

namespace Lizeral{
    class Mesh;
namespace Reflection{
namespace TypeFieldReflectionOparator{
    class TypeMeshOperator{
    public:
        static const char* getClassName(){ return "Mesh";}
        static void* constructorWithJson(const PJson& json_context){
            Mesh* ret_instance= new Mesh;
            PSerializer::read(json_context, *ret_instance);
            return ret_instance;
        }
        static PJson writeByName(void* instance){
            return PSerializer::write(*(Mesh*)instance);
        }
        // base class
        static int getMeshBaseClassReflectionInstanceList(ReflectionInstance* &out_list, void* instance){
            int count = 0;
            
            return count;
        }
        // fields
        
    };
}//namespace TypeFieldReflectionOparator


    void TypeWrapperRegister_Mesh(){
        
        
        
        ClassFunctionTuple* f_class_function_tuple_Mesh=new ClassFunctionTuple(
            &TypeFieldReflectionOparator::TypeMeshOperator::getMeshBaseClassReflectionInstanceList,
            &TypeFieldReflectionOparator::TypeMeshOperator::constructorWithJson,
            &TypeFieldReflectionOparator::TypeMeshOperator::writeByName);
        REGISTER_BASE_CLASS_TO_MAP("Mesh", f_class_function_tuple_Mesh);
    }
namespace TypeWrappersRegister{
        void Mesh(){ TypeWrapperRegister_Mesh();}
}//namespace TypeWrappersRegister

}//namespace Reflection
}//namespace Lizeral


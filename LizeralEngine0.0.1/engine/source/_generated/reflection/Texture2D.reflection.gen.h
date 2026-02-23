#pragma once
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
        static int getTexture2DBaseClassReflectionInstanceList(ReflectionInstance* &out_list, void* instance){
            int count = 1;
            out_list = new ReflectionInstance[count];
            for (int i=0;i<count;++i){
               out_list[i] = TypeMetaDef(Lizeral::Resource,static_cast<Texture2D*>(instance));
            }
            return count;
        }
        // fields
        
    };
}//namespace TypeFieldReflectionOparator


    void TypeWrapperRegister_Texture2D(){
        
        
        
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


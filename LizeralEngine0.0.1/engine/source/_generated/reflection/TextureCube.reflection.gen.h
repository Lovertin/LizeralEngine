#pragma once
#include "runtime\function\res_type\texture\TextureCube.h"

namespace Lizeral{
    class TextureCube;
namespace Reflection{
namespace TypeFieldReflectionOparator{
    class TypeTextureCubeOperator{
    public:
        static const char* getClassName(){ return "TextureCube";}
        static void* constructorWithJson(const PJson& json_context){
            TextureCube* ret_instance= new TextureCube;
            PSerializer::read(json_context, *ret_instance);
            return ret_instance;
        }
        static PJson writeByName(void* instance){
            return PSerializer::write(*(TextureCube*)instance);
        }
        // base class
        static int getTextureCubeBaseClassReflectionInstanceList(ReflectionInstance* &out_list, void* instance){
            int count = 1;
            out_list = new ReflectionInstance[count];
            for (int i=0;i<count;++i){
               out_list[i] = TypeMetaDef(Lizeral::Resource,static_cast<TextureCube*>(instance));
            }
            return count;
        }
        // fields
        
    };
}//namespace TypeFieldReflectionOparator


    void TypeWrapperRegister_TextureCube(){
        
        
        
        ClassFunctionTuple* f_class_function_tuple_TextureCube=new ClassFunctionTuple(
            &TypeFieldReflectionOparator::TypeTextureCubeOperator::getTextureCubeBaseClassReflectionInstanceList,
            &TypeFieldReflectionOparator::TypeTextureCubeOperator::constructorWithJson,
            &TypeFieldReflectionOparator::TypeTextureCubeOperator::writeByName);
        REGISTER_BASE_CLASS_TO_MAP("TextureCube", f_class_function_tuple_TextureCube);
    }
namespace TypeWrappersRegister{
        void TextureCube(){ TypeWrapperRegister_TextureCube();}
}//namespace TypeWrappersRegister

}//namespace Reflection
}//namespace Lizeral


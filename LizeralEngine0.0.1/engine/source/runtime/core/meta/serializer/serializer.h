#pragma once
#include "runtime/core/meta/json.h"
#include <cassert>
#include <type_traits>

// =====================================================================
// 1. 【关键修复】：前置声明 Reflection 命名空间中的类，彻底打破循环依赖！
// =====================================================================
namespace Lizeral {
    namespace Reflection {
        class TypeMeta;
        template<typename T> class ReflectionPtr;
    }
}

namespace Lizeral
{
    template<typename...>
    inline constexpr bool always_false = false;

    class PSerializer
    {
    public:
        // 遇到需要使用 Reflection 的函数，我们在这里【只写声明，不写实现】
        template<typename T> static PJson writePointer(T* instance);
        template<typename T> static T*& readPointer(const PJson& json_context, T*& instance);
        template<typename T> static PJson write(const Reflection::ReflectionPtr<T>& instance);
        template<typename T> static T*& read(const PJson& json_context, Reflection::ReflectionPtr<T>& instance);

        // 泛型 write (不需要 reflection.h，可以直接实现)
        template<typename T>
        static PJson write(const T& instance)
        {
            if constexpr (std::is_pointer<T>::value) {
                return writePointer((T)instance);
            } else if constexpr (std::is_enum_v<T>) {
                return PJson(static_cast<int>(instance));
            } else if constexpr (std::is_abstract_v<T>) {
                return PJson::object();
            } else {
                static_assert(always_false<T>, "PSerializer::write<T> has not been implemented yet!");
                return PJson();
            }
        }

        // 泛型 read (不需要 reflection.h，可以直接实现)
        template<typename T>
        static T& read(const PJson& json_context, T& instance)
        {
            if constexpr (std::is_pointer<T>::value) {
                return readPointer(json_context, instance);
            } else if constexpr (std::is_enum_v<T>) {
                instance = static_cast<T>(static_cast<int>(json_context.number_value()));
                return instance;
            } else if constexpr (std::is_abstract_v<T>) {
                return instance;
            } else {
                static_assert(always_false<T>, "PSerializer::read<T> has not been implemented yet!");
                return instance;
            }
        }
    }; // PSerializer 类的定义在这里安全结束！
} // namespace Lizeral 结束

// =====================================================================
// 2. 【关键修复】：等 PSerializer 完全定义好之后，再去包含 reflection.h
// 这样当 reflection.h 展开时，它就能完美认识 PSerializer 了！
// =====================================================================
#include "runtime/core/meta/reflection/reflection.h"

// =====================================================================
// 3. 在文件底部，补全之前那些需要用到 TypeMeta 的模板函数的实现
// =====================================================================
namespace Lizeral 
{
    template<typename T>
    inline PJson PSerializer::writePointer(T* instance) {
        return PJson::object {{"$typeName", PJson {"*"}}, {"$context", PSerializer::write(*instance)}};
    }

    template<typename T>
    inline T*& PSerializer::readPointer(const PJson& json_context, T*& instance) {
        assert(instance == nullptr);
        std::string type_name = json_context["$typeName"].string_value();
        assert(!type_name.empty());
        if ('*' == type_name[0]) {
            instance = new T;
            read(json_context["$context"], *instance);
        } else {
            instance = static_cast<T*>(
                Reflection::TypeMeta::newFromNameAndPJson(type_name, json_context["$context"]).m_instance);
        }
        return instance;
    }

    template<typename T>
    inline PJson PSerializer::write(const Reflection::ReflectionPtr<T>& instance) {
        T* instance_ptr = static_cast<T*>(instance.operator->());
        std::string type_name = instance.getTypeName();
        return PJson::object {{"$typeName", PJson(type_name)},
                              {"$context", Reflection::TypeMeta::writeByName(type_name, instance_ptr)}};
    }

    template<typename T>
    inline T*& PSerializer::read(const PJson& json_context, Reflection::ReflectionPtr<T>& instance) {
        std::string type_name = json_context["$typeName"].string_value();
        instance.setTypeName(type_name);
        return readPointer(json_context, instance.getPtrReference());
    }

    // 基础类型的特化声明
    template<> PJson PSerializer::write(const char& instance);
    template<> char& PSerializer::read(const PJson& json_context, char& instance);
    template<> PJson PSerializer::write(const int& instance);
    template<> int& PSerializer::read(const PJson& json_context, int& instance);
    template<> PJson PSerializer::write(const unsigned int& instance);
    template<> unsigned int& PSerializer::read(const PJson& json_context, unsigned int& instance);
    template<> PJson PSerializer::write(const float& instance);
    template<> float& PSerializer::read(const PJson& json_context, float& instance);
    template<> PJson PSerializer::write(const double& instance);
    template<> double& PSerializer::read(const PJson& json_context, double& instance);
    template<> PJson PSerializer::write(const bool& instance);
    template<> bool& PSerializer::read(const PJson& json_context, bool& instance);
    template<> PJson PSerializer::write(const std::string& instance);
    template<> std::string& PSerializer::read(const PJson& json_context, std::string& instance);

} // 完美结束！注意结尾不要有额外的 } 了
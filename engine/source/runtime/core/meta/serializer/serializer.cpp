#include "serializer.h"
#include "runtime/function/ecs/components/Collider/ColliderComponent.h"
#include <assert.h>
namespace Lizeral
{

    template<>
    PJson PSerializer::write(const char& instance)
    {
        return PJson(instance);
    }
    template<>
    char& PSerializer::read(const PJson& json_context, char& instance)
    {
        assert(json_context.is_number());
        return instance = json_context.number_value();
    }

    template<>
    PJson PSerializer::write(const int& instance)
    {
        return PJson(instance);
    }
    template<>
    int& PSerializer::read(const PJson& json_context, int& instance)
    {
        assert(json_context.is_number());
        return instance = static_cast<int>(json_context.number_value());
    }

    template<>
    PJson PSerializer::write(const unsigned int& instance)
    {
        return PJson(static_cast<int>(instance));
    }
    template<>
    unsigned int& PSerializer::read(const PJson& json_context, unsigned int& instance)
    {
        assert(json_context.is_number());
        return instance = static_cast<unsigned int>(json_context.number_value());
    }

    template<>
    PJson PSerializer::write(const float& instance)
    {
        return PJson(instance);
    }
    template<>
    float& PSerializer::read(const PJson& json_context, float& instance)
    {
        assert(json_context.is_number());
        return instance = static_cast<float>(json_context.number_value());
    }

    template<>
    PJson PSerializer::write(const double& instance)
    {
        return PJson(instance);
    }
    template<>
    double& PSerializer::read(const PJson& json_context, double& instance)
    {
        assert(json_context.is_number());
        return instance = static_cast<float>(json_context.number_value());
    }

    template<>
    PJson PSerializer::write(const bool& instance)
    {
        return PJson(instance);
    }
    template<>
    bool& PSerializer::read(const PJson& json_context, bool& instance)
    {
        assert(json_context.is_bool());
        return instance = json_context.bool_value();
    }

    template<>
    PJson PSerializer::write(const std::string& instance)
    {
        return PJson(instance);
    }
    template<>
    std::string& PSerializer::read(const PJson& json_context, std::string& instance)
    {
        assert(json_context.is_string());
        return instance = json_context.string_value();
    }


    // // 2. 解决 Resource (抽象基类) 的序列化报错
    // template<>
    // inline json11::Json PSerializer::write(const Resource& instance) {
    //     return json11::Json::object(); // 返回空 JSON 对象
    // }
    // template<>
    // inline Resource& PSerializer::read(const json11::Json& json_context, Resource& instance) {
    //     return instance; // 啥也不干，直接返回
    // }

    // // 3. 解决 Material 的序列化报错
    // template<>
    // inline json11::Json PSerializer::write(const Material& instance) {
    //     return json11::Json::object();
    // }
    // template<>
    // inline Material& PSerializer::read(const json11::Json& json_context, Material& instance) {
    //     return instance;
    // }

    // template<>
    // PJson PSerializer::write(const Reflection::object& instance)
    //{
    // return PJson::object();
    //}
    // template<>
    // Reflection::object& PSerializer::read(const PJson& json_context, Reflection::object& instance)
    //{
    //	return instance;
    //}

    //////////////////////////////////
    // template of generation coder
    //////////////////////////////////

    /*template<>
    PJson PSerializer::write(const ss& instance)
    {
        return PJson();
    }
    template<>
    PJson PSerializer::write(const jkj& instance)
    {
        return PJson();
    }

    template<>
    PJson PSerializer::write(const test_class& instance)
    {
        PJson::array aa;
        for(auto& item: instance.c)
        {
            aa.emplace_back(PSerializer::write(item));
        }
        ss jj;
        reflection::object* jjj1 = &jj;
        auto kkk = PSerializer::write(jjj1);
        auto jjj = kkk.dump();
        return PJson::object{
            {"a",PSerializer::write(instance.a)},
            {"b",PSerializer::write(instance.b)},
            {"c",aa}
        };
    }
    template<>
    test_class& PSerializer::read(const PJson& json_context, test_class& instance)
    {
        assert(json_context.is_object());
        PSerializer::read(json_context["a"], instance.a);
        PSerializer::read(json_context["b"], instance.b);

        assert(json_context["c"].is_array());
        PJson::array cc = json_context["c"].array_items();
        instance.c.resize(cc.size());
        for (size_t index=0; index < cc.size();++index)
        {
            PSerializer::read(cc[index], instance.c[index]);
        }
        return instance;
    }*/

    ////////////////////////////////////
} // namespace Lizeral
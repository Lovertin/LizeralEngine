#pragma once
#include "ICommand.h"
#include "runtime/core/meta/reflection/reflection.h"
#include "runtime/function/ecs/entity.h"
#include <string>
#include <functional>

namespace Lizeral {

    template<typename T>
    class EditPropertyCommand : public ICommand {
    public:
        // 定义一个能动态获取组件安全地址的函数指针类型
        using ComponentGetter = std::function<void*()>;

        EditPropertyCommand(Entity entity, 
                            const std::string& compName, 
                            const std::string& fieldName, 
                            T oldValue, 
                            T newValue,
                            ComponentGetter getter)
            : m_Entity(entity), m_CompName(compName), m_FieldName(fieldName), 
              m_OldValue(oldValue), m_NewValue(newValue), m_CompGetter(getter) {}

        void Execute() override { ApplyValue(m_NewValue); }
        void Undo() override    { ApplyValue(m_OldValue); }

        bool MergeWith(const ICommand* other) override {
            const auto* otherCmd = dynamic_cast<const EditPropertyCommand<T>*>(other);
            if (!otherCmd) return false;

            // 必须是同一个实体、同一个组件的同一个字段
            if (m_Entity == otherCmd->m_Entity && 
                m_CompName == otherCmd->m_CompName && 
                m_FieldName == otherCmd->m_FieldName) 
            {
                // 保留旧值，吸收最新的目标值
                this->m_NewValue = otherCmd->m_NewValue;
                return true;
            }
            return false;
        }

    private:
        void ApplyValue(T& value) {
            // 1. 极其关键：每次都动态去 ECS 查询内存地址！绝对不会野指针！
            void* compInstance = m_CompGetter(); 
            if (!compInstance) return;

            // 2. 利用你的反射系统拿到访问器
            Reflection::TypeMeta meta = Reflection::TypeMeta::newMetaFromName(m_CompName);
            if (meta.isValid()) {
                Reflection::FieldAccessor field = meta.getFieldByName(m_FieldName.c_str());
                
                // 3. 写入数据。由于你的反射生成代码中已经带有 TryNotifyReflectionUpdated，
                // 这里不仅会改值，还会完美触发你的 DIRTY 标记！
                field.set(compInstance, &value); 
            }
        }

    private:
        Entity m_Entity;
        std::string m_CompName;
        std::string m_FieldName;
        T m_OldValue;
        T m_NewValue;
        ComponentGetter m_CompGetter; // 安全获取指针的闭包
    };

} // namespace Lizeral
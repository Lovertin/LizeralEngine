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

            if (m_Entity == otherCmd->m_Entity && 
                m_CompName == otherCmd->m_CompName && 
                m_FieldName == otherCmd->m_FieldName) 
            {
                this->m_NewValue = otherCmd->m_NewValue;
                return true;
            }
            return false;
        }

    private:
        void ApplyValue(T& value) {
            void* compInstance = m_CompGetter(); 
            if (!compInstance) return;

            Reflection::TypeMeta meta = Reflection::TypeMeta::newMetaFromName(m_CompName);
            if (meta.isValid()) {
                Reflection::FieldAccessor field = meta.getFieldByName(m_FieldName.c_str());

                field.set(compInstance, &value); 
            }
        }

    private:
        Entity m_Entity;
        std::string m_CompName;
        std::string m_FieldName;
        T m_OldValue;
        T m_NewValue;
        ComponentGetter m_CompGetter;
    };

} // namespace Lizeral4
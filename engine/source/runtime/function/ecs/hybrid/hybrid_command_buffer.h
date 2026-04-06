#pragma once
#include "runtime/function/ecs/hybrid/hybrid_registry.h"
#include <memory>
#include <utility>
#include <vector>

namespace Lizeral {

    class HybridCommandBuffer {
    private:
        struct ICommand {
            virtual ~ICommand() = default;
            virtual void apply(HybridRegistry& registry) = 0;
        };

        template<typename T>
        struct EmplaceCommand final : ICommand {
            Entity entity;
            T component;

            EmplaceCommand(Entity target, T value)
                : entity(target), component(std::move(value)) {}

            void apply(HybridRegistry& registry) override {
                registry.emplace<T>(entity, std::move(component));
            }
        };

        template<typename T>
        struct RemoveCommand final : ICommand {
            Entity entity;

            explicit RemoveCommand(Entity target)
                : entity(target) {}

            void apply(HybridRegistry& registry) override {
                registry.remove<T>(entity);
            }
        };

        struct DestroyCommand final : ICommand {
            Entity entity;

            explicit DestroyCommand(Entity target)
                : entity(target) {}

            void apply(HybridRegistry& registry) override {
                registry.destroy(entity);
            }
        };

    public:
        template<typename T>
        void emplace(Entity entity, T component) {
            m_commands.push_back(std::make_unique<EmplaceCommand<T>>(entity, std::move(component)));
        }

        template<typename T>
        void remove(Entity entity) {
            m_commands.push_back(std::make_unique<RemoveCommand<T>>(entity));
        }

        void destroy(Entity entity) {
            m_commands.push_back(std::make_unique<DestroyCommand>(entity));
        }

        void flush(HybridRegistry& registry) {
            for (const std::unique_ptr<ICommand>& command : m_commands) {
                command->apply(registry);
            }

            m_commands.clear();
        }

        void clear() {
            m_commands.clear();
        }

        bool empty() const {
            return m_commands.empty();
        }

    private:
        std::vector<std::unique_ptr<ICommand>> m_commands;
    };

} // namespace Lizeral

#pragma once
#include "entity.h"
#include <vector>
#include <unordered_map>
#include <cassert>
#include <iostream>

namespace Lizeral {

    class IPool {
    public:
        virtual ~IPool() = default;
        virtual bool has(Entity entity) const = 0;
        virtual void remove(Entity entity) = 0;

        virtual void clear() = 0;
    };

    // Sparse Set
    template<typename T>
    class Pool : public IPool {
    private:
        std::vector<T>      m_components; 
        
        // Swap & Pop
        std::vector<Entity> m_dense_entities; 

        // Entity ID -> index 
        std::unordered_map<Entity, size_t> m_entity_to_index;

    public:

        bool has(Entity entity) const override {
            return m_entity_to_index.find(entity) != m_entity_to_index.end();
        }

        // emplace component
        template<typename... Args>
        T& emplace(Entity entity, Args&&... args) {
            assert(!has(entity));

            m_components.emplace_back(std::forward<Args>(args)...);
            m_dense_entities.push_back(entity);

            m_entity_to_index[entity] = m_components.size() - 1;

            return m_components.back();
        }

        T& get(Entity entity) {
            assert(has(entity));
            size_t index = m_entity_to_index[entity];
            return m_components[index];
        }

        void remove(Entity entity) override {
            if (!has(entity)) return;

            size_t removed_index = m_entity_to_index[entity];
            size_t last_index = m_components.size() - 1;
            Entity last_entity = m_dense_entities.back();

            if (removed_index != last_index) {
                m_components[removed_index] = std::move(m_components.back());
                m_dense_entities[removed_index] = last_entity;

                m_entity_to_index[last_entity] = removed_index;
            }

            m_components.pop_back();
            m_dense_entities.pop_back();
            m_entity_to_index.erase(entity);
        }

        void clear() override {
            m_components.clear();
            m_dense_entities.clear();
            m_entity_to_index.clear();
        }

        auto begin() { return m_components.begin(); }
        auto end() { return m_components.end(); }

        const std::vector<Entity>& get_entities() const { return m_dense_entities; }
    };
}
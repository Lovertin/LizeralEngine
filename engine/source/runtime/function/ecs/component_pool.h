#pragma once
#include "entity.h"
#include <limits>
#include <vector>
#include <cassert>

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
        static constexpr size_t npos = std::numeric_limits<size_t>::max();

        std::vector<T>      m_components; 
        
        // Swap & Pop
        std::vector<Entity> m_dense_entities; 

        // Entity ID -> dense index
        std::vector<size_t> m_sparse_indices;

    public:

        bool has(Entity entity) const override {
            const size_t entity_index = static_cast<size_t>(entity);
            if (entity_index >= m_sparse_indices.size()) {
                return false;
            }

            const size_t dense_index = m_sparse_indices[entity_index];
            return dense_index != npos &&
                   dense_index < m_dense_entities.size() &&
                   m_dense_entities[dense_index] == entity;
        }

        // emplace component
        template<typename... Args>
        T& emplace(Entity entity, Args&&... args) {
            assert(!has(entity));

            const size_t entity_index = static_cast<size_t>(entity);
            if (entity_index >= m_sparse_indices.size()) {
                m_sparse_indices.resize(entity_index + 1, npos);
            }

            m_components.emplace_back(std::forward<Args>(args)...);
            m_dense_entities.push_back(entity);

            m_sparse_indices[entity_index] = m_components.size() - 1;

            return m_components.back();
        }

        T& get(Entity entity) {
            assert(has(entity));
            return m_components[m_sparse_indices[static_cast<size_t>(entity)]];
        }

        const T& get(Entity entity) const {
            assert(has(entity));
            return m_components[m_sparse_indices[static_cast<size_t>(entity)]];
        }

        void remove(Entity entity) override {
            if (!has(entity)) return;

            size_t removed_index = m_sparse_indices[static_cast<size_t>(entity)];
            size_t last_index = m_components.size() - 1;
            Entity last_entity = m_dense_entities.back();

            if (removed_index != last_index) {
                m_components[removed_index] = std::move(m_components.back());
                m_dense_entities[removed_index] = last_entity;

                m_sparse_indices[static_cast<size_t>(last_entity)] = removed_index;
            }

            m_components.pop_back();
            m_dense_entities.pop_back();
            m_sparse_indices[static_cast<size_t>(entity)] = npos;
        }

        void clear() override {
            m_components.clear();
            m_dense_entities.clear();
            m_sparse_indices.clear();
        }

        auto begin() { return m_components.begin(); }
        auto end() { return m_components.end(); }
        auto begin() const { return m_components.begin(); }
        auto end() const { return m_components.end(); }

        size_t size() const { return m_components.size(); }
        const std::vector<Entity>& get_entities() const { return m_dense_entities; }
    };
}

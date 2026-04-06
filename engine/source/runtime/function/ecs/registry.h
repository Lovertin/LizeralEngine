#pragma once
#include "entity.h"
#include "component_pool.h"
#include "view.h"
#include <unordered_map>
#include <typeindex>
#include <memory>
#include <algorithm>
#include <limits>
#include <tuple>

namespace Lizeral {

    class Registry {
    private:
        Entity m_next_entity_id = 0;
        
        std::unordered_map<std::type_index, std::shared_ptr<IPool>> m_pools;

    public:
        Entity create() {
            return m_next_entity_id++;
        }

        template<typename T>
        Pool<T>* get_pool() {
            std::type_index type = std::type_index(typeid(T));
            if (m_pools.find(type) == m_pools.end()) {
                m_pools[type] = std::make_shared<Pool<T>>();
            }
            return static_cast<Pool<T>*>(m_pools[type].get());
        }

        template<typename T>
        Pool<T>* try_get_pool() {
            std::type_index type = std::type_index(typeid(T));
            auto it = m_pools.find(type);
            if (it == m_pools.end()) {
                return nullptr;
            }

            return static_cast<Pool<T>*>(it->second.get());
        }

        template<typename T>
        const Pool<T>* try_get_pool() const {
            std::type_index type = std::type_index(typeid(T));
            auto it = m_pools.find(type);
            if (it == m_pools.end()) {
                return nullptr;
            }

            return static_cast<const Pool<T>*>(it->second.get());
        }

        template<typename T, typename... Args>
        T& emplace(Entity entity, Args&&... args) {
            return get_pool<T>()->emplace(entity, std::forward<Args>(args)...);
        }

        template<typename T>
        T& get(Entity entity) {
            auto* pool = try_get_pool<T>();
            assert(pool && "Component pool does not exist.");
            return pool->get(entity);
        }

        template<typename T>
        const T& get(Entity entity) const {
            auto* pool = try_get_pool<T>();
            assert(pool && "Component pool does not exist.");
            return pool->get(entity);
        }

        template<typename T>
        bool has(Entity entity) {
            auto* pool = try_get_pool<T>();
            return pool && pool->has(entity);
        }

        template<typename T>
        bool has(Entity entity) const {
            auto* pool = try_get_pool<T>();
            return pool && pool->has(entity);
        }

        template<typename T>
        void remove(Entity entity) {
            auto* pool = try_get_pool<T>();
            if (pool) {
                pool->remove(entity);
            }
        }

        void destroy(Entity entity) {
            for (auto& pair : m_pools) {
                pair.second->remove(entity);
            }
        }

        void clear() {

            for (auto& pair : m_pools) {
                if (pair.second) {
                    pair.second->clear();
                }
            }

            m_next_entity_id = 0; 
        }

        template<typename... Components>
        View<Components...> view() {
            static const std::vector<Entity> empty_entities;

            const std::vector<Entity>* smallest_pool_entities = &empty_entities;
            size_t smallest_pool_size = std::numeric_limits<size_t>::max();
            bool has_missing_pool = false;

            ([&] {
                auto* pool = try_get_pool<Components>();
                if (!pool) {
                    has_missing_pool = true;
                    return;
                }

                if (pool->size() < smallest_pool_size) {
                    smallest_pool_size = pool->size();
                    smallest_pool_entities = &pool->get_entities();
                }
            }(), ...);

            if (has_missing_pool) {
                return View<Components...>(this, &empty_entities);
            }

            return View<Components...>(this, smallest_pool_entities);
        }
    };
    
    template<typename... Components>
    void View<Components...>::Iterator::advance_until_valid() {

        while (index < entities.size()) { 
            Entity current = entities[index];
            // Fold expression (C++17)
            if ((registry->has<Components>(current) && ...)) {
                return; // Found valid entity
            }
            index++;
        }
    }

    template<typename... Components>
    typename View<Components...>::Iterator View<Components...>::begin() {
        return Iterator(m_registry, *m_smallest_set, 0);
    }

    template<typename... Components>
    typename View<Components...>::Iterator View<Components...>::end() {
        return Iterator(m_registry, *m_smallest_set, m_smallest_set->size());
    }

    template<typename... Components>
    template<typename T>
    T& View<Components...>::get(Entity entity) {
        return m_registry->get<T>(entity);
    }
}

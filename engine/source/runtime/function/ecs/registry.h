#pragma once
#include "entity.h"
#include "component_pool.h"
#include "view.h"
#include <unordered_map>
#include <typeindex>
#include <memory>
#include <algorithm>

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

        template<typename T, typename... Args>
        T& emplace(Entity entity, Args&&... args) {
            return get_pool<T>()->emplace(entity, std::forward<Args>(args)...);
        }

        template<typename T>
        T& get(Entity entity) {
            return get_pool<T>()->get(entity);
        }

        template<typename T>
        bool has(Entity entity) {
            return get_pool<T>()->has(entity);
        }

        template<typename T>
        void remove(Entity entity) {
            get_pool<T>()->remove(entity);
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
            auto* pool = get_pool<std::tuple_element_t<0, std::tuple<Components...>>>();
            return View<Components...>(this, &pool->get_entities());
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
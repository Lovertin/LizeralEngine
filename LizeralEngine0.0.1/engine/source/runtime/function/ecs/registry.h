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
        
        // 存储所有类型的组件池
        // Key: 组件类型的 type_index, Value: Pool<T>
        std::unordered_map<std::type_index, std::shared_ptr<IPool>> m_pools;

    public:
        Entity create() {
            return m_next_entity_id++;
        }

        // 获取或创建某个类型的池子
        template<typename T>
        Pool<T>* get_pool() {
            std::type_index type = std::type_index(typeid(T));
            if (m_pools.find(type) == m_pools.end()) {
                m_pools[type] = std::make_shared<Pool<T>>();
            }
            return static_cast<Pool<T>*>(m_pools[type].get());
        }

        // 添加组件
        template<typename T, typename... Args>
        T& emplace(Entity entity, Args&&... args) {
            return get_pool<T>()->emplace(entity, std::forward<Args>(args)...);
        }

        // 获取组件
        template<typename T>
        T& get(Entity entity) {
            return get_pool<T>()->get(entity);
        }
        
        // 检查是否有组件
        template<typename T>
        bool has(Entity entity) {
            return get_pool<T>()->has(entity);
        }

        template<typename T>
        void remove(Entity entity) {
            get_pool<T>()->remove(entity);
        }

        // 【新增】：彻底销毁实体！遍历所有组件池并将其拔除
        void destroy(Entity entity) {
            for (auto& pair : m_pools) {
                pair.second->remove(entity);
            }
        }

        // 核心：创建视图
        template<typename... Components>
        View<Components...> view() {
            // 简单优化：找到实体数量最少的那个组件池作为主  理，直接拿第一个类型的池子遍历
            // (更优的做法是比较 sizeof m_dense_entities 并选最小的)
            auto* pool = get_pool<std::tuple_element_t<0, std::tuple<Components...>>>();
            return View<Components...>(this, &pool->get_entities());
        }
    };

    // --- View 的模板实现 ---
    
    template<typename... Components>
    void View<Components...>::Iterator::advance_until_valid() {
        // 遍历 entities，直到找到一个 entity 同时拥有所有 Components
        while (index < entities.size()) { 
            Entity current = entities[index];
            // Fold expression (C++17): 检查 registry 是否拥有所有组件
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
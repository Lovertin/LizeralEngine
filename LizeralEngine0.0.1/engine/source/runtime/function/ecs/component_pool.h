#pragma once
#include "entity.h"
#include <vector>
#include <unordered_map>
#include <cassert>
#include <iostream>

namespace Lizeral {

    // 1. 定义一个虚基类，方便 Registry 统一管理不同类型的池子
    class IPool {
    public:
        virtual ~IPool() = default;
        virtual bool has(Entity entity) const = 0;
        virtual void remove(Entity entity) = 0;
    };

    // 2. 具体的组件池 (Sparse Set 实现)
    template<typename T>
    class Pool : public IPool {
    private:
        // 紧密排列的数据，CPU 缓存非常友好
        std::vector<T>      m_components; 
        
        // 紧密排列的 Entity ID，与 m_components 一一对应
        // 用于从 index 反查 Entity，或者在删除时做 Swap & Pop
        std::vector<Entity> m_dense_entities; 

        // 稀疏映射: Entity ID -> index (在 m_components 中的下标)
        // 也可以用 std::vector (如果 ID 较小) 或 std::unordered_map (如果 ID 很稀疏)
        std::unordered_map<Entity, size_t> m_entity_to_index;

    public:

        bool has(Entity entity) const override {
            return m_entity_to_index.find(entity) != m_entity_to_index.end();
        }

        // 添加组件
        template<typename... Args>
        T& emplace(Entity entity, Args&&... args) {
            assert(!has(entity));

            // 1. 在 Dense 尾部构造组件
            m_components.emplace_back(std::forward<Args>(args)...);
            m_dense_entities.push_back(entity);

            // 2. 记录映射关系
            m_entity_to_index[entity] = m_components.size() - 1;

            return m_components.back();
        }

        // 获取组件
        T& get(Entity entity) {
            assert(has(entity));
            size_t index = m_entity_to_index[entity];
            return m_components[index];
        }

        // 删除组件 (Swap & Pop 技巧：O(1) 删除)
        void remove(Entity entity) override {
            if (!has(entity)) return;

            size_t removed_index = m_entity_to_index[entity];
            size_t last_index = m_components.size() - 1;
            Entity last_entity = m_dense_entities.back();

            // 1. 将最后一个元素移动到被删除的位置
            if (removed_index != last_index) {
                m_components[removed_index] = std::move(m_components.back());
                m_dense_entities[removed_index] = last_entity;
                
                // 更新最后一个元素的索引指向
                m_entity_to_index[last_entity] = removed_index;
            }

            // 2. 移除末尾
            m_components.pop_back();
            m_dense_entities.pop_back();
            m_entity_to_index.erase(entity);
        }

        // 提供给 View 遍历用的迭代器
        auto begin() { return m_components.begin(); }
        auto end() { return m_components.end(); }
        
        // 获取所有拥有该组件的 Entity 列表 (用于 View 逻辑)
        const std::vector<Entity>& get_entities() const { return m_dense_entities; }
    };
}
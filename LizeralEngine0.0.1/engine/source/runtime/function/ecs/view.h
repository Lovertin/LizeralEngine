#pragma once
#include "component_pool.h"

namespace Lizeral {

    // 简单的 View，用于遍历单个或多个组件
    template<typename... Components>
    class View {
    private:
        // 保存所有相关组件池的指针
        // 实际上我们只需要遍历数量最少的那个池子，然后去其他池子检查 has()
        class Registry* m_registry;
        const std::vector<Entity>* m_smallest_set;

    public:
        View(Registry* registry, const std::vector<Entity>* smallest) 
            : m_registry(registry), m_smallest_set(smallest) {}

        // 自定义迭代器，用于支持 range-based for loop
        struct Iterator {
            Registry* registry;
            const std::vector<Entity>& entities;
            size_t index;

            Iterator(Registry* r, const std::vector<Entity>& e, size_t i) 
                : registry(r), entities(e), index(i) {
                // 初始化时通过有效性检查找到第一个合法的
                advance_until_valid();
            }

            void advance_until_valid(); // 在 Registry.cpp 或下方实现

            Entity operator*() const { return entities[index]; }
            
            Iterator& operator++() {
                index++;
                advance_until_valid();
                return *this;
            }

            bool operator!=(const Iterator& other) const { return index != other.index; }
        };

        Iterator begin();
        Iterator end();
        
        // 方便获取组件的辅助函数: view.get<T>(entity)
        template<typename T>
        T& get(Entity entity);
    };
}
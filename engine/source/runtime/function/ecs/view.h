#pragma once
#include "component_pool.h"

namespace Lizeral {

    template<typename... Components>
    class View {
    private:

        class Registry* m_registry;
        const std::vector<Entity>* m_smallest_set;

    public:
        View(Registry* registry, const std::vector<Entity>* smallest) 
            : m_registry(registry), m_smallest_set(smallest) {}

        struct Iterator {
            Registry* registry;
            const std::vector<Entity>& entities;
            size_t index;

            Iterator(Registry* r, const std::vector<Entity>& e, size_t i) 
                : registry(r), entities(e), index(i) {
                advance_until_valid();
            }

            void advance_until_valid(); 

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

        template<typename T>
        T& get(Entity entity);
    };
}
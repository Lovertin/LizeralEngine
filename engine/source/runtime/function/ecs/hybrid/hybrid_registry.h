#pragma once
#include "runtime/function/ecs/entity.h"
#include "runtime/function/ecs/component_pool.h"
#include "runtime/function/ecs/hybrid/hybrid_component_traits.h"
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <new>
#include <type_traits>
#include <typeindex>
#include <unordered_map>
#include <utility>
#include <vector>

namespace Lizeral {

    struct HybridRegistryConfig {
        size_t chunk_capacity = 64;
    };

    class HybridRegistry {
    public:
        using ComponentTypeId = std::uint32_t;
        using ArchetypeSignature = std::vector<ComponentTypeId>;

        struct ComponentTypeInfo {
            ComponentTypeId id = 0;
            size_t size = 0;
            size_t alignment = 0;
            void (*move_construct)(void* dst, void* src) = nullptr;
            void (*destroy)(void* instance) = nullptr;
        };

        template<typename... Components>
        class ChunkView;

        template<typename... Components>
        class View;

    private:
        struct SignatureHash {
            size_t operator()(const ArchetypeSignature& signature) const noexcept {
                size_t hash = 1469598103934665603ull;
                for (ComponentTypeId id : signature) {
                    hash ^= static_cast<size_t>(id) + 0x9e3779b97f4a7c15ull + (hash << 6) + (hash >> 2);
                }
                return hash;
            }
        };

        struct ChunkColumn {
            const ComponentTypeInfo* info = nullptr;
            void* data = nullptr;

            ChunkColumn() = default;

            ChunkColumn(const ComponentTypeInfo& type_info, size_t capacity)
                : info(&type_info) {
                data = ::operator new(type_info.size * capacity, std::align_val_t(type_info.alignment));
            }

            ChunkColumn(const ChunkColumn&) = delete;
            ChunkColumn& operator=(const ChunkColumn&) = delete;

            ChunkColumn(ChunkColumn&& other) noexcept
                : info(other.info), data(other.data) {
                other.info = nullptr;
                other.data = nullptr;
            }

            ChunkColumn& operator=(ChunkColumn&& other) noexcept {
                if (this == &other) {
                    return *this;
                }

                release();
                info = other.info;
                data = other.data;
                other.info = nullptr;
                other.data = nullptr;
                return *this;
            }

            ~ChunkColumn() {
                release();
            }

            void* slot(size_t row) {
                return static_cast<std::byte*>(data) + (row * info->size);
            }

            const void* slot(size_t row) const {
                return static_cast<const std::byte*>(data) + (row * info->size);
            }

        private:
            void release() {
                if (data && info) {
                    ::operator delete(data, std::align_val_t(info->alignment));
                }
                data = nullptr;
                info = nullptr;
            }
        };

        class Archetype;

        class Chunk {
        public:
            Chunk(Archetype* owner, const std::vector<const ComponentTypeInfo*>& component_infos, size_t capacity)
                : m_owner(owner), m_capacity(capacity), m_entities(capacity, null_entity) {
                m_columns.reserve(component_infos.size());
                for (const ComponentTypeInfo* info : component_infos) {
                    m_columns.emplace_back(*info, capacity);
                    m_column_by_type.emplace(info->id, m_columns.size() - 1);
                }
            }

            Chunk(const Chunk&) = delete;
            Chunk& operator=(const Chunk&) = delete;

            ~Chunk() {
                destroy_all();
            }

            bool has_space() const { return m_size < m_capacity; }
            size_t size() const { return m_size; }
            size_t capacity() const { return m_capacity; }

            size_t allocate(Entity entity) {
                assert(has_space());
                const size_t row = m_size++;
                m_entities[row] = entity;
                return row;
            }

            Entity entity_at(size_t row) const {
                assert(row < m_size);
                return m_entities[row];
            }

            const Entity* entity_data() const {
                return m_entities.data();
            }

            bool has_component(ComponentTypeId type_id) const {
                return m_column_by_type.find(type_id) != m_column_by_type.end();
            }

            void* get_component_ptr(ComponentTypeId type_id, size_t row) {
                assert(row < m_size || row < m_capacity);
                return m_columns[m_column_by_type.at(type_id)].slot(row);
            }

            const void* get_component_ptr(ComponentTypeId type_id, size_t row) const {
                assert(row < m_size);
                return m_columns.at(m_column_by_type.at(type_id)).slot(row);
            }

            template<typename T>
            T& get(size_t row) {
                return *static_cast<T*>(get_component_ptr(get_component_type_info<T>().id, row));
            }

            template<typename T>
            const T& get(size_t row) const {
                return *static_cast<const T*>(get_component_ptr(get_component_type_info<T>().id, row));
            }

            template<typename MovedEntityCallback>
            void erase_row(size_t row, MovedEntityCallback&& moved_callback) {
                assert(row < m_size);

                const size_t last_row = m_size - 1;
                if (row != last_row) {
                    const Entity moved_entity = m_entities[last_row];
                    for (ChunkColumn& column : m_columns) {
                        column.info->destroy(column.slot(row));
                        column.info->move_construct(column.slot(row), column.slot(last_row));
                        column.info->destroy(column.slot(last_row));
                    }
                    m_entities[row] = moved_entity;
                    moved_callback(moved_entity, row);
                } else {
                    for (ChunkColumn& column : m_columns) {
                        column.info->destroy(column.slot(row));
                    }
                }

                m_entities[last_row] = null_entity;
                --m_size;
            }

        private:
            void destroy_all() {
                for (size_t row = 0; row < m_size; ++row) {
                    for (ChunkColumn& column : m_columns) {
                        column.info->destroy(column.slot(row));
                    }
                }
                m_size = 0;
            }

        private:
            Archetype* m_owner = nullptr;
            size_t m_capacity = 0;
            size_t m_size = 0;
            std::vector<Entity> m_entities;
            std::vector<ChunkColumn> m_columns;
            std::unordered_map<ComponentTypeId, size_t> m_column_by_type;
        };

        class Archetype {
        public:
            Archetype(ArchetypeSignature signature, std::vector<const ComponentTypeInfo*> component_infos, size_t chunk_capacity)
                : m_signature(std::move(signature)),
                  m_component_infos(std::move(component_infos)),
                  m_chunk_capacity(chunk_capacity) {
                for (size_t i = 0; i < m_component_infos.size(); ++i) {
                    m_component_index.emplace(m_component_infos[i]->id, i);
                }
            }

            const ArchetypeSignature& signature() const { return m_signature; }
            const std::vector<const ComponentTypeInfo*>& component_infos() const { return m_component_infos; }

            bool contains(ComponentTypeId type_id) const {
                return m_component_index.find(type_id) != m_component_index.end();
            }

            bool matches(const ArchetypeSignature& query_signature) const {
                return std::includes(
                    m_signature.begin(), m_signature.end(),
                    query_signature.begin(), query_signature.end());
            }

            Chunk& get_or_create_chunk_with_space() {
                for (const std::unique_ptr<Chunk>& chunk : m_chunks) {
                    if (chunk->has_space()) {
                        return *chunk;
                    }
                }

                m_chunks.push_back(std::make_unique<Chunk>(this, m_component_infos, m_chunk_capacity));
                return *m_chunks.back();
            }

            const std::vector<std::unique_ptr<Chunk>>& chunks() const {
                return m_chunks;
            }

        private:
            ArchetypeSignature m_signature;
            std::vector<const ComponentTypeInfo*> m_component_infos;
            std::unordered_map<ComponentTypeId, size_t> m_component_index;
            std::vector<std::unique_ptr<Chunk>> m_chunks;
            size_t m_chunk_capacity = 0;
        };

        struct EntityRecord {
            bool alive = false;
            Archetype* archetype = nullptr;
            Chunk* chunk = nullptr;
            size_t row = 0;
        };

        struct QueryCacheEntry {
            size_t structure_version = 0;
            std::vector<Archetype*> archetypes;
        };

    public:
        explicit HybridRegistry(HybridRegistryConfig config = {})
            : m_config(config) {
            assert(m_config.chunk_capacity > 0);
        }

        Entity create() {
            const Entity entity = m_next_entity_id++;
            if (static_cast<size_t>(entity) >= m_entity_records.size()) {
                m_entity_records.resize(static_cast<size_t>(entity) + 1);
            }

            EntityRecord& record = m_entity_records[static_cast<size_t>(entity)];
            record = {};
            record.alive = true;
            return entity;
        }

        bool is_alive(Entity entity) const {
            const size_t index = static_cast<size_t>(entity);
            return index < m_entity_records.size() && m_entity_records[index].alive;
        }

        template<typename T, typename... Args>
        T& emplace(Entity entity, Args&&... args) {
            require_entity(entity);

            if constexpr (is_hybrid_archetype_component_v<T>) {
                return emplace_archetype<T>(entity, std::forward<Args>(args)...);
            } else {
                return get_or_create_sparse_pool<T>()->emplace(entity, std::forward<Args>(args)...);
            }
        }

        template<typename T>
        bool has(Entity entity) const {
            if (!is_alive(entity)) {
                return false;
            }

            if constexpr (is_hybrid_archetype_component_v<T>) {
                const EntityRecord& record = m_entity_records[static_cast<size_t>(entity)];
                return record.archetype != nullptr &&
                       record.chunk != nullptr &&
                       record.archetype->contains(get_component_type_info<T>().id);
            } else {
                const Pool<T>* pool = try_get_sparse_pool<T>();
                return pool && pool->has(entity);
            }
        }

        template<typename T>
        T& get(Entity entity) {
            T* component = try_get<T>(entity);
            assert(component && "Component does not exist on entity.");
            return *component;
        }

        template<typename T>
        const T& get(Entity entity) const {
            const T* component = try_get<T>(entity);
            assert(component && "Component does not exist on entity.");
            return *component;
        }

        template<typename T>
        T* try_get(Entity entity) {
            if (!has<T>(entity)) {
                return nullptr;
            }

            if constexpr (is_hybrid_archetype_component_v<T>) {
                EntityRecord& record = m_entity_records[static_cast<size_t>(entity)];
                return &record.chunk->template get<T>(record.row);
            } else {
                Pool<T>* pool = try_get_sparse_pool<T>();
                return pool ? &pool->get(entity) : nullptr;
            }
        }

        template<typename T>
        const T* try_get(Entity entity) const {
            if (!has<T>(entity)) {
                return nullptr;
            }

            if constexpr (is_hybrid_archetype_component_v<T>) {
                const EntityRecord& record = m_entity_records[static_cast<size_t>(entity)];
                return &record.chunk->template get<T>(record.row);
            } else {
                const Pool<T>* pool = try_get_sparse_pool<T>();
                return pool ? &pool->get(entity) : nullptr;
            }
        }

        template<typename T>
        void remove(Entity entity) {
            if (!is_alive(entity) || !has<T>(entity)) {
                return;
            }

            if constexpr (is_hybrid_archetype_component_v<T>) {
                remove_archetype<T>(entity);
            } else {
                Pool<T>* pool = try_get_sparse_pool<T>();
                if (pool) {
                    pool->remove(entity);
                }
            }
        }

        void destroy(Entity entity) {
            if (!is_alive(entity)) {
                return;
            }

            EntityRecord& record = m_entity_records[static_cast<size_t>(entity)];
            if (record.chunk) {
                Chunk* source_chunk = record.chunk;
                const size_t source_row = record.row;
                source_chunk->erase_row(source_row, [this, source_chunk](Entity moved_entity, size_t new_row) {
                    EntityRecord& moved_record = m_entity_records[static_cast<size_t>(moved_entity)];
                    moved_record.chunk = source_chunk;
                    moved_record.row = new_row;
                });
            }

            for (auto& pair : m_sparse_pools) {
                pair.second->remove(entity);
            }

            record = {};
        }

        void clear() {
            m_sparse_pools.clear();
            m_archetypes.clear();
            m_query_cache.clear();
            m_entity_records.clear();
            m_next_entity_id = 0;
            ++m_structure_version;
        }

        template<typename... Components>
        View<Components...> view() {
            if constexpr ((is_hybrid_archetype_component_v<Components> && ...)) {
                return View<Components...>(this, collect_matching_chunks<Components...>());
            } else {
                return View<Components...>(this, collect_matching_entities<Components...>());
            }
        }

        template<typename... Components>
        std::vector<ChunkView<Components...>> chunk_query() {
            static_assert((is_hybrid_archetype_component_v<Components> && ...),
                          "Chunk queries are only supported for archetype-managed components.");

            std::vector<ChunkView<Components...>> views;
            std::vector<Chunk*> chunks = collect_matching_chunks<Components...>();
            views.reserve(chunks.size());
            for (Chunk* chunk : chunks) {
                views.emplace_back(this, chunk);
            }

            return views;
        }

        template<typename... Components, typename Func>
        void for_each_chunk(Func&& func) {
            static_assert((is_hybrid_archetype_component_v<Components> && ...),
                          "Chunk iteration is only supported for archetype-managed components.");

            for (ChunkView<Components...>& chunk : chunk_query<Components...>()) {
                func(chunk);
            }
        }

    private:
        template<typename T>
        static const ComponentTypeInfo& get_component_type_info() {
            static_assert(std::is_move_constructible_v<T>,
                          "Archetype-managed components must be move constructible.");
            static_assert(std::is_destructible_v<T>,
                          "Archetype-managed components must be destructible.");

            static const ComponentTypeInfo info = {
                next_component_type_id(),
                sizeof(T),
                alignof(T),
                [](void* dst, void* src) {
                    new (dst) T(std::move(*static_cast<T*>(src)));
                },
                [](void* instance) {
                    static_cast<T*>(instance)->~T();
                }
            };

            return info;
        }

        static ComponentTypeId next_component_type_id() {
            static ComponentTypeId next_id = 0;
            return next_id++;
        }

        EntityRecord& require_entity(Entity entity) {
            assert(is_alive(entity) && "Entity is not alive in HybridRegistry.");
            return m_entity_records[static_cast<size_t>(entity)];
        }

        template<typename T>
        Pool<T>* get_or_create_sparse_pool() {
            const std::type_index type = std::type_index(typeid(T));
            auto it = m_sparse_pools.find(type);
            if (it == m_sparse_pools.end()) {
                it = m_sparse_pools.emplace(type, std::make_shared<Pool<T>>()).first;
            }

            return static_cast<Pool<T>*>(it->second.get());
        }

        template<typename T>
        Pool<T>* try_get_sparse_pool() {
            const std::type_index type = std::type_index(typeid(T));
            auto it = m_sparse_pools.find(type);
            return it == m_sparse_pools.end() ? nullptr : static_cast<Pool<T>*>(it->second.get());
        }

        template<typename T>
        const Pool<T>* try_get_sparse_pool() const {
            const std::type_index type = std::type_index(typeid(T));
            auto it = m_sparse_pools.find(type);
            return it == m_sparse_pools.end() ? nullptr : static_cast<const Pool<T>*>(it->second.get());
        }

        static ArchetypeSignature add_component_to_signature(ArchetypeSignature signature, ComponentTypeId type_id) {
            auto it = std::lower_bound(signature.begin(), signature.end(), type_id);
            signature.insert(it, type_id);
            return signature;
        }

        static ArchetypeSignature remove_component_from_signature(ArchetypeSignature signature, ComponentTypeId type_id) {
            auto it = std::lower_bound(signature.begin(), signature.end(), type_id);
            if (it != signature.end() && *it == type_id) {
                signature.erase(it);
            }
            return signature;
        }

        Archetype* get_or_create_archetype(const ArchetypeSignature& signature) {
            auto it = m_archetypes.find(signature);
            if (it != m_archetypes.end()) {
                return it->second.get();
            }

            std::vector<const ComponentTypeInfo*> component_infos;
            component_infos.reserve(signature.size());
            for (ComponentTypeId id : signature) {
                component_infos.push_back(m_component_info_by_id.at(id));
            }

            auto archetype = std::make_unique<Archetype>(signature, std::move(component_infos), m_config.chunk_capacity);
            Archetype* archetype_ptr = archetype.get();
            m_archetypes.emplace(signature, std::move(archetype));
            ++m_structure_version;
            return archetype_ptr;
        }

        template<typename T>
        void register_archetype_component_type() {
            const ComponentTypeInfo& info = get_component_type_info<T>();
            m_component_info_by_id.emplace(info.id, &info);
        }

        template<typename T, typename... Args>
        T& emplace_archetype(Entity entity, Args&&... args) {
            register_archetype_component_type<T>();

            EntityRecord& record = require_entity(entity);
            assert(!has<T>(entity) && "Component already exists on entity.");

            const ComponentTypeId new_type_id = get_component_type_info<T>().id;
            const ArchetypeSignature source_signature =
                record.archetype ? record.archetype->signature() : ArchetypeSignature{};
            const ArchetypeSignature target_signature = add_component_to_signature(source_signature, new_type_id);

            Archetype* target_archetype = get_or_create_archetype(target_signature);
            Chunk& target_chunk = target_archetype->get_or_create_chunk_with_space();
            const size_t target_row = target_chunk.allocate(entity);

            if (record.chunk) {
                for (const ComponentTypeInfo* info : target_archetype->component_infos()) {
                    void* dst = target_chunk.get_component_ptr(info->id, target_row);
                    if (info->id == new_type_id) {
                        new (dst) T(std::forward<Args>(args)...);
                    } else {
                        info->move_construct(dst, record.chunk->get_component_ptr(info->id, record.row));
                    }
                }
            } else {
                void* dst = target_chunk.get_component_ptr(new_type_id, target_row);
                new (dst) T(std::forward<Args>(args)...);
            }

            Chunk* source_chunk = record.chunk;
            const size_t source_row = record.row;
            record.archetype = target_archetype;
            record.chunk = &target_chunk;
            record.row = target_row;

            if (source_chunk) {
                source_chunk->erase_row(source_row, [this, source_chunk](Entity moved_entity, size_t new_row) {
                    EntityRecord& moved_record = m_entity_records[static_cast<size_t>(moved_entity)];
                    moved_record.chunk = source_chunk;
                    moved_record.row = new_row;
                });
            }

            return target_chunk.template get<T>(target_row);
        }

        template<typename T>
        void remove_archetype(Entity entity) {
            EntityRecord& record = require_entity(entity);
            assert(record.archetype && record.chunk);

            const ComponentTypeId remove_type_id = get_component_type_info<T>().id;
            const ArchetypeSignature target_signature =
                remove_component_from_signature(record.archetype->signature(), remove_type_id);

            Archetype* target_archetype = nullptr;
            Chunk* target_chunk = nullptr;
            size_t target_row = 0;

            if (!target_signature.empty()) {
                target_archetype = get_or_create_archetype(target_signature);
                target_chunk = &target_archetype->get_or_create_chunk_with_space();
                target_row = target_chunk->allocate(entity);

                for (const ComponentTypeInfo* info : target_archetype->component_infos()) {
                    void* dst = target_chunk->get_component_ptr(info->id, target_row);
                    info->move_construct(dst, record.chunk->get_component_ptr(info->id, record.row));
                }
            }

            Chunk* source_chunk = record.chunk;
            const size_t source_row = record.row;
            record.archetype = target_archetype;
            record.chunk = target_chunk;
            record.row = target_row;

            source_chunk->erase_row(source_row, [this, source_chunk](Entity moved_entity, size_t new_row) {
                EntityRecord& moved_record = m_entity_records[static_cast<size_t>(moved_entity)];
                moved_record.chunk = source_chunk;
                moved_record.row = new_row;
            });
        }

        template<typename... Components>
        ArchetypeSignature build_query_signature() {
            ArchetypeSignature signature;
            (((is_hybrid_archetype_component_v<Components>)
                  ? (register_archetype_component_type<Components>(),
                     signature.push_back(get_component_type_info<Components>().id), 0)
                  : 0), ...);
            std::sort(signature.begin(), signature.end());
            return signature;
        }

        template<typename... Components>
        std::vector<Chunk*> collect_matching_chunks() {
            std::vector<Chunk*> chunks;
            const ArchetypeSignature query_signature = build_query_signature<Components...>();
            for (Archetype* archetype : get_cached_matching_archetypes(query_signature)) {
                for (const std::unique_ptr<Chunk>& chunk : archetype->chunks()) {
                    if (chunk->size() > 0) {
                        chunks.push_back(chunk.get());
                    }
                }
            }

            return chunks;
        }

        const std::vector<Archetype*>& get_cached_matching_archetypes(const ArchetypeSignature& query_signature) {
            auto cache_it = m_query_cache.find(query_signature);
            if (cache_it != m_query_cache.end() && cache_it->second.structure_version == m_structure_version) {
                return cache_it->second.archetypes;
            }

            QueryCacheEntry entry{};
            entry.structure_version = m_structure_version;
            for (auto& pair : m_archetypes) {
                Archetype* archetype = pair.second.get();
                if (archetype->matches(query_signature)) {
                    entry.archetypes.push_back(archetype);
                }
            }

            auto [it, inserted] = m_query_cache.insert_or_assign(query_signature, std::move(entry));
            (void)inserted;
            return it->second.archetypes;
        }

        template<typename T>
        bool update_sparse_driver(const std::vector<Entity>*& best_entities, size_t& best_size) {
            if constexpr (!is_hybrid_archetype_component_v<T>) {
                const Pool<T>* pool = try_get_sparse_pool<T>();
                if (!pool) {
                    return false;
                }

                if (pool->size() < best_size) {
                    best_size = pool->size();
                    best_entities = &pool->get_entities();
                }
            }

            return true;
        }

        template<typename... Components>
        std::vector<Entity> collect_matching_entities() {
            std::vector<Entity> entities;

            const std::vector<Entity>* driver_entities = nullptr;
            size_t driver_size = static_cast<size_t>(-1);
            const bool has_driver = (update_sparse_driver<Components>(driver_entities, driver_size) && ...);
            if (!has_driver || !driver_entities) {
                return entities;
            }

            entities.reserve(driver_entities->size());
            for (Entity entity : *driver_entities) {
                if ((has<Components>(entity) && ...)) {
                    entities.push_back(entity);
                }
            }

            return entities;
        }

    private:
        HybridRegistryConfig m_config;
        Entity m_next_entity_id = 0;
        std::vector<EntityRecord> m_entity_records;
        std::unordered_map<std::type_index, std::shared_ptr<IPool>> m_sparse_pools;
        std::unordered_map<ComponentTypeId, const ComponentTypeInfo*> m_component_info_by_id;
        std::unordered_map<ArchetypeSignature, std::unique_ptr<Archetype>, SignatureHash> m_archetypes;
        std::unordered_map<ArchetypeSignature, QueryCacheEntry, SignatureHash> m_query_cache;
        size_t m_structure_version = 1;

    public:
        template<typename... Components>
        class ChunkView {
        public:
            ChunkView(HybridRegistry* registry, Chunk* chunk)
                : m_registry(registry), m_chunk(chunk) {}

            size_t size() const { return m_chunk ? m_chunk->size() : 0; }

            const Entity* entities() const {
                return m_chunk ? m_chunk->entity_data() : nullptr;
            }

            template<typename T>
            T& get(size_t row) {
                return m_chunk->template get<T>(row);
            }

        private:
            HybridRegistry* m_registry = nullptr;
            Chunk* m_chunk = nullptr;
        };

        template<typename... Components>
        class View {
        public:
            View(HybridRegistry* registry, std::vector<Chunk*> chunks)
                : m_registry(registry), m_chunks(std::move(chunks)), m_chunk_mode(true) {}

            View(HybridRegistry* registry, std::vector<Entity> entities)
                : m_registry(registry), m_entities(std::move(entities)), m_chunk_mode(false) {}

            struct Iterator {
                View* owner = nullptr;
                size_t outer_index = 0;
                size_t inner_index = 0;

                Iterator(View* view, size_t outer, size_t inner)
                    : owner(view), outer_index(outer), inner_index(inner) {
                    advance_to_valid();
                }

                void advance_to_valid() {
                    if (!owner || !owner->m_chunk_mode) {
                        return;
                    }

                    while (outer_index < owner->m_chunks.size()) {
                        Chunk* chunk = owner->m_chunks[outer_index];
                        if (inner_index < chunk->size()) {
                            return;
                        }

                        ++outer_index;
                        inner_index = 0;
                    }
                }

                Entity operator*() const {
                    if (owner->m_chunk_mode) {
                        return owner->m_chunks[outer_index]->entity_at(inner_index);
                    }

                    return owner->m_entities[outer_index];
                }

                Iterator& operator++() {
                    if (owner->m_chunk_mode) {
                        ++inner_index;
                        advance_to_valid();
                    } else {
                        ++outer_index;
                    }
                    return *this;
                }

                bool operator!=(const Iterator& other) const {
                    return outer_index != other.outer_index || inner_index != other.inner_index || owner != other.owner;
                }
            };

            Iterator begin() {
                return Iterator(this, 0, 0);
            }

            Iterator end() {
                return m_chunk_mode ? Iterator(this, m_chunks.size(), 0) : Iterator(this, m_entities.size(), 0);
            }

            template<typename T>
            T& get(Entity entity) {
                return m_registry->template get<T>(entity);
            }

        private:
            HybridRegistry* m_registry = nullptr;
            std::vector<Chunk*> m_chunks;
            std::vector<Entity> m_entities;
            bool m_chunk_mode = false;
        };
    };

} // namespace Lizeral

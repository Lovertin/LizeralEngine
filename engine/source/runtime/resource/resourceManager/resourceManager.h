// LizeralEngine0.0.1\engine\source\runtime\resourceManager\resourceManager.h
#pragma once
#include <unordered_map>
#include <string>
#include <memory>
#include <functional>
#include <mutex>
#include <fstream>
#include "runtime/resource/resource.h"

namespace Lizeral {
    // 资源加载优先级
    enum class ResourcePriority {
        LOW,       
        NORMAL,    
        HIGH,     
        CRITICAL   
    };
    
    using ResourceLoadedCallback = std::function<void(std::shared_ptr<ResourceAsset>)>;
    
    class ResourceManager {
    public:
        static ResourceManager& GetInstance() {
            static ResourceManager instance;
            return instance;
        }

        template<typename T>
        std::shared_ptr<T> Load(const std::string& path, 
                               ResourcePriority priority = ResourcePriority::NORMAL) {

            std::string fullPath = ResolvePath(path);
            {
                std::lock_guard<std::mutex> lock(m_cacheMutex);
                auto it = m_cache.find(path);
                if (it != m_cache.end()) {
                    return std::static_pointer_cast<T>(it->second);
                }
            }
            std::shared_ptr<T> resource = std::make_shared<T>();
            resource->SetPath(path);
            
            if (resource->LoadFromFile(path)) {
                std::lock_guard<std::mutex> lock(m_cacheMutex);
                m_cache[path] = resource;
                
                auto callbackIt = m_callbacks.find(path);
                if (callbackIt != m_callbacks.end()) {
                    callbackIt->second(resource);
                    m_callbacks.erase(callbackIt);
                }
                
                return resource;
            }

            return nullptr;
        }
        
        template<typename T>
        void LoadAsync(const std::string& path, 
                      ResourceLoadedCallback callback = nullptr,
                      ResourcePriority priority = ResourcePriority::NORMAL) {
            std::string fullPath = ResolvePath(path);
            if (callback) {
                std::lock_guard<std::mutex> lock(m_cacheMutex);
                m_callbacks[fullPath] = callback;
            }
            
            auto resource = Load<T>(path, priority);
            // if (callback && resource) {
            //     callback(resource);
            // }
        }
 
        template<typename T>
        bool Preload(const std::string& path) {
            return Load<T>(path, ResourcePriority::LOW) != nullptr;
        }
        
        void Unload(const std::string& path) {
            std::lock_guard<std::mutex> lock(m_cacheMutex);
            m_cache.erase(path);
        }

        void UnloadUnusedResources() {
            std::lock_guard<std::mutex> lock(m_cacheMutex);
            for (auto it = m_cache.begin(); it != m_cache.end(); ) {

                if (it->second.use_count() == 1) {
                    it = m_cache.erase(it);
                } else {
                    ++it;
                }
            }
        }
        
        void ClearAll() {
            std::lock_guard<std::mutex> lock(m_cacheMutex);
            m_cache.clear();
            m_callbacks.clear();
        }
        bool IsLoaded(const std::string& path) const {
            std::lock_guard<std::mutex> lock(m_cacheMutex);
            return m_cache.find(path) != m_cache.end();
        }
        
        size_t GetResourceSize(const std::string& path) const {
            return 0;
        }
        
        size_t GetTotalMemoryUsage() const {
            return 0;
        }
        void AddSearchPath(const std::string& path) {
            m_searchPaths.push_back(path);
        }
        
        void SetRootPath(const std::string& path) {
            m_rootPath = path;
        }

    private:
        ResourceManager() = default;
        ~ResourceManager() = default;
        
        ResourceManager(const ResourceManager&) = delete;
        ResourceManager& operator=(const ResourceManager&) = delete;
        
        std::string ResolvePath(const std::string& path) const {
            if (path.empty() || path[0] == '/' || path[0] == '\\' || 
                (path.size() > 1 && path[1] == ':')) {
                return path;
            }
            
            for (const auto& searchPath : m_searchPaths) {
                std::string fullPath = m_rootPath + "/" + searchPath + "/" + path;
                std::ifstream f(fullPath);
                if (f.good()) {
                    return fullPath;
                }
            }
            
            return m_rootPath + "/" + path;
        }
        
    private:
        mutable std::mutex m_cacheMutex;
        std::unordered_map<std::string, std::shared_ptr<ResourceAsset>> m_cache;
        std::unordered_map<std::string, ResourceLoadedCallback> m_callbacks;
        
        std::string m_rootPath = "assets";  
        std::vector<std::string> m_searchPaths = {
            "models", "textures", "shaders", "materials"
        };
    };
}

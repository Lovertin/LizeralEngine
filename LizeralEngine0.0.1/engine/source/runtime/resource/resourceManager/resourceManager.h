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
        LOW,        // 后台加载
        NORMAL,     // 普通优先级
        HIGH,       // 立即加载
        CRITICAL    // 阻塞加载
    };
    
    // 资源加载回调
    using ResourceLoadedCallback = std::function<void(std::shared_ptr<Resource>)>;
    
    class ResourceManager {
    public:
        static ResourceManager& GetInstance() {
            static ResourceManager instance;
            return instance;
        }

        // 基础加载接口（同步）
        template<typename T>
        std::shared_ptr<T> Load(const std::string& path, 
                               ResourcePriority priority = ResourcePriority::NORMAL) {

            //将标准路径转化成绝对路径
            std::string fullPath = ResolvePath(path);
            // 1. 查缓存（线程安全）
            {
                std::lock_guard<std::mutex> lock(m_cacheMutex);
                auto it = m_cache.find(path);
                if (it != m_cache.end()) {
                    return std::static_pointer_cast<T>(it->second);
                }
            }

            // 2. 缓存未命中，创建新实例
            std::shared_ptr<T> resource = std::make_shared<T>();
            resource->SetPath(path);
            
            // 3. 加载资源
            if (resource->LoadFromFile(path)) {
                std::lock_guard<std::mutex> lock(m_cacheMutex);
                m_cache[path] = resource;
                
                // 触发加载完成回调
                auto callbackIt = m_callbacks.find(path);
                if (callbackIt != m_callbacks.end()) {
                    callbackIt->second(resource);
                    m_callbacks.erase(callbackIt);
                }
                
                return resource;
            }

            return nullptr;
        }
        
        // 异步加载接口
        template<typename T>
        void LoadAsync(const std::string& path, 
                      ResourceLoadedCallback callback = nullptr,
                      ResourcePriority priority = ResourcePriority::NORMAL) {
            std::string fullPath = ResolvePath(path);
            if (callback) {
                std::lock_guard<std::mutex> lock(m_cacheMutex);
                m_callbacks[fullPath] = callback;
            }
            
            // 这里可以启动一个后台线程来加载
            // 简化版本直接同步加载
            auto resource = Load<T>(path, priority);
            // if (callback && resource) {
            //     callback(resource);
            // }
        }
        
        // 预加载资源（不返回，只加入缓存）
        template<typename T>
        bool Preload(const std::string& path) {
            return Load<T>(path, ResourcePriority::LOW) != nullptr;
        }
        
        // 卸载资源
        void Unload(const std::string& path) {
            std::lock_guard<std::mutex> lock(m_cacheMutex);
            m_cache.erase(path);
        }

        void UnloadUnusedResources() {
            std::lock_guard<std::mutex> lock(m_cacheMutex);
            for (auto it = m_cache.begin(); it != m_cache.end(); ) {
                // use_count() == 1 说明除了 ResourceManager 的 m_cache 字典拿着它，
                // 没有任何 Component 或外界变量还在用它。此时可以安全卸载显存！
                if (it->second.use_count() == 1) {
                    it = m_cache.erase(it);
                } else {
                    ++it;
                }
            }
        }
        
        // 清空所有资源
        void ClearAll() {
            std::lock_guard<std::mutex> lock(m_cacheMutex);
            m_cache.clear();
            m_callbacks.clear();
        }
        
        // 检查资源是否已加载
        bool IsLoaded(const std::string& path) const {
            std::lock_guard<std::mutex> lock(m_cacheMutex);
            return m_cache.find(path) != m_cache.end();
        }
        
        // 获取资源大小（内存占用）
        size_t GetResourceSize(const std::string& path) const {
            // 需要资源类提供获取大小的接口
            return 0;
        }
        
        // 获取总内存占用
        size_t GetTotalMemoryUsage() const {
            // 遍历所有资源计算大小
            return 0;
        }
        
        // 添加资源搜索路径
        void AddSearchPath(const std::string& path) {
            m_searchPaths.push_back(path);
        }
        
        // 设置资源根目录
        void SetRootPath(const std::string& path) {
            m_rootPath = path;
        }

    private:
        ResourceManager() = default;
        ~ResourceManager() = default;
        
        // 禁止拷贝
        ResourceManager(const ResourceManager&) = delete;
        ResourceManager& operator=(const ResourceManager&) = delete;
        
        // 辅助函数：解析完整路径
        std::string ResolvePath(const std::string& path) const {
            // 如果已经是绝对路径，直接返回
            if (path.empty() || path[0] == '/' || path[0] == '\\' || 
                (path.size() > 1 && path[1] == ':')) {
                return path;
            }
            
            // 在搜索路径中查找文件
            for (const auto& searchPath : m_searchPaths) {
                std::string fullPath = m_rootPath + "/" + searchPath + "/" + path;
                // 检查文件是否存在
                std::ifstream f(fullPath);
                if (f.good()) {
                    return fullPath;
                }
            }
            
            // 默认返回根路径下的文件
            return m_rootPath + "/" + path;
        }
        
    private:
        mutable std::mutex m_cacheMutex;
        std::unordered_map<std::string, std::shared_ptr<Resource>> m_cache;
        std::unordered_map<std::string, ResourceLoadedCallback> m_callbacks;
        
        std::string m_rootPath = "assets";  // 资源根目录
        std::vector<std::string> m_searchPaths = {
            "models", "textures", "shaders", "materials"
        };
    };
}
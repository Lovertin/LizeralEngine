// source/runtime/function/render/rhi/RHIContext.h
#pragma once
#include <string>
#include <vector>

namespace Lizeral {

    class IRHIContext {
    public:
        virtual ~IRHIContext() = default;

        // 初始化图形 API (参数是你的应用名称)
        virtual void Initialize(const std::string& appName,const std::vector<const char*>& windowExtensions = {}) = 0;
        
        // 销毁并清理资源
        virtual void Shutdown() = 0;

        // 获取底层原始指针 (只有 RHI 内部的强转才用得到)
        virtual void* GetNativeInstance() const = 0; 
    };

} // namespace Lizeral
// source/runtime/function/render/rhi/RHIContext.h
#pragma once
#include <string>
#include <vector>

namespace Lizeral {

    class IRHIContext {
    public:
        virtual ~IRHIContext() = default;

        virtual void Initialize(const std::string& appName,const std::vector<const char*>& windowExtensions = {}) = 0;
        
        virtual void Shutdown() = 0;

        virtual void* GetNativeInstance() const = 0; 
    };

} // namespace Lizeral
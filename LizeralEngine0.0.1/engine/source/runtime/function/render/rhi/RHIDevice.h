// source/runtime/function/render/rhi/RHIDevice.h
#pragma once

namespace Lizeral {

    class IRHIDevice {
    public:
        virtual ~IRHIDevice() = default;

        virtual void WaitIdle() = 0;
    };

} // namespace Lizeral
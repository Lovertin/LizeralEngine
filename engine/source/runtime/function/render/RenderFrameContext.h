#pragma once

#include "runtime/core/math/matrix4.h"
#include "runtime/core/math/vector3.h"
#include <vulkan/vulkan.h>
#include <cstdint>

namespace Lizeral {

    struct RenderFrameContext {
        VkCommandBuffer cmd { VK_NULL_HANDLE };
        VkViewport viewport {};
        VkRect2D scissor {};
        uint32_t width { 0 };
        uint32_t height { 0 };
        uint32_t frameIndex { 0 };

        Matrix4x4 view;
        Matrix4x4 projection;
        Matrix4x4 viewProj;

        Vector3 cameraPos;
        Vector3 cameraRight { 1.0f, 0.0f, 0.0f };
        Vector3 cameraUp { 0.0f, 1.0f, 0.0f };
        Vector3 cameraForward { 0.0f, 0.0f, -1.0f };

        VkDescriptorSet sceneDepthSet { VK_NULL_HANDLE };
    };

} // namespace Lizeral

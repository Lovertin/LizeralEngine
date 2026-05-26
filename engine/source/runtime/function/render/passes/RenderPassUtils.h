#pragma once

#include <vulkan/vulkan.h>

namespace Lizeral {

    void TransitionImageLayout(
        VkCommandBuffer cmd,
        VkImage image,
        VkImageLayout oldLayout,
        VkImageLayout newLayout,
        VkImageAspectFlags aspectMask
    );

} // namespace Lizeral

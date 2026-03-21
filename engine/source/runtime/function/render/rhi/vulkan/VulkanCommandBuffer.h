#pragma once
#include <vulkan/vulkan.h>

namespace Lizeral {

    class VulkanDevice;
    class VulkanCommandPool;

    class VulkanCommandBuffer {
    public:
        VulkanCommandBuffer(VulkanDevice* device, VulkanCommandPool* pool);
        ~VulkanCommandBuffer();

        void Begin(VkCommandBufferUsageFlags flags = 0);
        void End();
        
        void SubmitAndIdle();

        VkCommandBuffer GetNativeBuffer() const { return m_commandBuffer; }

    private:
        VulkanDevice* m_device { nullptr };
        VulkanCommandPool* m_pool { nullptr };
        VkCommandBuffer m_commandBuffer { VK_NULL_HANDLE };
    };

} // namespace Lizeral
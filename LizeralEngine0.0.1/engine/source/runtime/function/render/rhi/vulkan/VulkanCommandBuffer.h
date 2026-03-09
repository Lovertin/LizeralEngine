#pragma once
#include <vulkan/vulkan.h>

namespace Lizeral {

    class VulkanDevice;
    class VulkanCommandPool;

    class VulkanCommandBuffer {
    public:
        // 构造时直接从 Pool 里分配出一张清单
        VulkanCommandBuffer(VulkanDevice* device, VulkanCommandPool* pool);
        ~VulkanCommandBuffer();

        // --- 核心录制流程 ---
        void Begin(VkCommandBufferUsageFlags flags = 0);
        void End();
        
        // 发车！并阻塞 CPU 等待 GPU 把这张清单上的活儿干完
        void SubmitAndIdle();

        VkCommandBuffer GetNativeBuffer() const { return m_commandBuffer; }

    private:
        VulkanDevice* m_device { nullptr };
        VulkanCommandPool* m_pool { nullptr };
        VkCommandBuffer m_commandBuffer { VK_NULL_HANDLE };
    };

} // namespace Lizeral
#pragma once
#include <vulkan/vulkan.h>
#include <string>
#include <vector>

namespace Lizeral {

    class VulkanDevice;

    class VulkanComputePipeline {
    public:
        VulkanComputePipeline(VulkanDevice* device, const std::string& shaderPath);
        ~VulkanComputePipeline();

        VulkanComputePipeline(const VulkanComputePipeline&) = delete;
        VulkanComputePipeline& operator=(const VulkanComputePipeline&) = delete;

        VkPipeline GetNativePipeline() const { return m_pipeline; }
        VkPipelineLayout GetPipelineLayout() const { return m_pipelineLayout; }

    private:
        VulkanDevice* m_device { nullptr };
        VkPipeline m_pipeline { VK_NULL_HANDLE };
        VkPipelineLayout m_pipelineLayout { VK_NULL_HANDLE };

        std::vector<char> readFile(const std::string& filename);
        VkShaderModule createShaderModule(const std::vector<char>& code);
    };

} // namespace Lizeral
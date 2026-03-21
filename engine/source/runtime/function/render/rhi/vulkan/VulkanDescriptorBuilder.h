#pragma once

#include <vulkan/vulkan.h>
#include <vector>
#include <unordered_map>
#include <deque>

namespace Lizeral {

    class VulkanDevice;

    class VulkanDescriptorBuilder {
    public:
        VulkanDescriptorBuilder();

        VulkanDescriptorBuilder& BindImage(uint32_t binding, const VkDescriptorImageInfo* imageInfo, VkDescriptorType type, VkShaderStageFlags stageFlags);
        
        VulkanDescriptorBuilder& BindImageArray(uint32_t binding, const VkDescriptorImageInfo* imageInfos, uint32_t layoutCount, uint32_t writeCount, VkDescriptorType type, VkShaderStageFlags stageFlags, bool isBindless = false);

        VulkanDescriptorBuilder& BindBuffer(uint32_t binding, const VkDescriptorBufferInfo* bufferInfo, VkDescriptorType type, VkShaderStageFlags stageFlags);

        VulkanDescriptorBuilder& BindAccelerationStructure(uint32_t binding, const VkAccelerationStructureKHR* as, VkShaderStageFlags stageFlags);

        bool Build(VulkanDevice* device, VkDescriptorSetLayout& outLayout, VkDescriptorPool& outPool, VkDescriptorSet& outSet);

        bool Build(VulkanDevice* device, VkDescriptorPool pool, VkDescriptorSetLayout& outLayout, VkDescriptorSet& outSet);

    private:
        std::vector<VkDescriptorSetLayoutBinding> m_bindings;
        std::vector<VkDescriptorBindingFlags> m_bindingFlags;
        std::vector<VkWriteDescriptorSet> m_writes;
        
        std::deque<VkWriteDescriptorSetAccelerationStructureKHR> m_asWrites;

        bool m_useBindless;
    };

} // namespace Lizeral
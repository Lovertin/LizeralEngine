#pragma once

#include <vulkan/vulkan.h>
#include <vector>
#include <unordered_map>

namespace Lizeral {

    class VulkanDevice;

    class VulkanDescriptorBuilder {
    public:
        VulkanDescriptorBuilder();

        // 1. 绑定单张贴图 (常规用途)
        VulkanDescriptorBuilder& BindImage(uint32_t binding, const VkDescriptorImageInfo* imageInfo, VkDescriptorType type, VkShaderStageFlags stageFlags);
        
        // 2. 绑定贴图数组 (专为 Bindless 打造！)
        VulkanDescriptorBuilder& BindImageArray(uint32_t binding, const VkDescriptorImageInfo* imageInfos, uint32_t layoutCount, uint32_t writeCount, VkDescriptorType type, VkShaderStageFlags stageFlags, bool isBindless = false);

        // 3. 绑定 Buffer (如果将来需要传 UBO / SSBO 可以用到)
        VulkanDescriptorBuilder& BindBuffer(uint32_t binding, const VkDescriptorBufferInfo* bufferInfo, VkDescriptorType type, VkShaderStageFlags stageFlags);

        // ★ 核心构建器 A：包含 Layout、分配专属 Pool，并直接完成数据的 Write
        bool Build(VulkanDevice* device, VkDescriptorSetLayout& outLayout, VkDescriptorPool& outPool, VkDescriptorSet& outSet);

        // ★ 核心构建器 B：复用已有的 Pool
        bool Build(VulkanDevice* device, VkDescriptorPool pool, VkDescriptorSetLayout& outLayout, VkDescriptorSet& outSet);

    private:
        std::vector<VkDescriptorSetLayoutBinding> m_bindings;
        std::vector<VkDescriptorBindingFlags> m_bindingFlags;
        std::vector<VkWriteDescriptorSet> m_writes;
        
        bool m_useBindless;
    };

} // namespace Lizeral
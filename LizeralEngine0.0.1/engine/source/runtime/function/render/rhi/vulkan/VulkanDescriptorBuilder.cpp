#include "VulkanDescriptorBuilder.h"
#include "VulkanDevice.h"
#include <stdexcept>
#include <iostream>

namespace Lizeral {

    VulkanDescriptorBuilder::VulkanDescriptorBuilder() : m_useBindless(false) {}

    VulkanDescriptorBuilder& VulkanDescriptorBuilder::BindImageArray(
        uint32_t binding, const VkDescriptorImageInfo* imageInfos, 
        uint32_t layoutCount, uint32_t writeCount, 
        VkDescriptorType type, VkShaderStageFlags stageFlags, bool isBindless) 
    {
        // 1. 记录布局绑定
        VkDescriptorSetLayoutBinding layoutBinding{};
        layoutBinding.binding = binding;
        layoutBinding.descriptorCount = layoutCount;
        layoutBinding.descriptorType = type;
        layoutBinding.stageFlags = stageFlags;
        m_bindings.push_back(layoutBinding);

        // 2. 记录绑定标志 (如果是 Bindless，允许未填满)
        if (isBindless) {
            m_bindingFlags.push_back(VK_DESCRIPTOR_BINDING_PARTIALLY_BOUND_BIT);
            m_useBindless = true;
        } else {
            m_bindingFlags.push_back(0); // 普通情况不加标志
        }

        // 3. 记录写入操作 (利用指针，推迟到 Build() 时再执行)
        VkWriteDescriptorSet write{VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET};
        write.dstBinding = binding;
        write.dstArrayElement = 0;
        write.descriptorType = type;
        write.descriptorCount = writeCount; // 实际写入的数量
        write.pImageInfo = imageInfos;
        m_writes.push_back(write);

        return *this;
    }

    VulkanDescriptorBuilder& VulkanDescriptorBuilder::BindImage(
        uint32_t binding, const VkDescriptorImageInfo* imageInfo, 
        VkDescriptorType type, VkShaderStageFlags stageFlags) 
    {
        // 单张图片就是特殊的数组，数量为 1
        return BindImageArray(binding, imageInfo, 1, 1, type, stageFlags, false);
    }

    VulkanDescriptorBuilder& VulkanDescriptorBuilder::BindBuffer(
        uint32_t binding, const VkDescriptorBufferInfo* bufferInfo, 
        VkDescriptorType type, VkShaderStageFlags stageFlags) 
    {
        VkDescriptorSetLayoutBinding layoutBinding{};
        layoutBinding.binding = binding;
        layoutBinding.descriptorCount = 1;
        layoutBinding.descriptorType = type;
        layoutBinding.stageFlags = stageFlags;
        m_bindings.push_back(layoutBinding);

        m_bindingFlags.push_back(0);

        VkWriteDescriptorSet write{VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET};
        write.dstBinding = binding;
        write.dstArrayElement = 0;
        write.descriptorType = type;
        write.descriptorCount = 1;
        write.pBufferInfo = bufferInfo;
        m_writes.push_back(write);

        return *this;
    }

    bool VulkanDescriptorBuilder::Build(VulkanDevice* device, VkDescriptorSetLayout& outLayout, VkDescriptorPool& outPool, VkDescriptorSet& outSet) {
        // 1. 自动统计当前需要的所有 PoolSize
        std::unordered_map<VkDescriptorType, uint32_t> poolSizeMap;
        for (const auto& binding : m_bindings) {
            poolSizeMap[binding.descriptorType] += binding.descriptorCount;
        }

        std::vector<VkDescriptorPoolSize> poolSizes;
        for (const auto& pair : poolSizeMap) {
            poolSizes.push_back({ pair.first, pair.second });
        }

        // 2. 创建专属池子
        VkDescriptorPoolCreateInfo poolInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO};
        poolInfo.poolSizeCount = static_cast<uint32_t>(poolSizes.size());
        poolInfo.pPoolSizes = poolSizes.data();
        poolInfo.maxSets = 1; // 这个 Builder 一次只建一个 Set

        if (vkCreateDescriptorPool(device->GetNativeDevice(), &poolInfo, nullptr, &outPool) != VK_SUCCESS) {
            throw std::runtime_error("Failed to create descriptor pool via Builder!");
        }

        // 3. 复用另一个 Build 函数完成剩余工作
        return Build(device, outPool, outLayout, outSet);
    }

    bool VulkanDescriptorBuilder::Build(VulkanDevice* device, VkDescriptorPool pool, VkDescriptorSetLayout& outLayout, VkDescriptorSet& outSet) {
        // 1. 创建 Layout (根据是否为 Bindless 决定是否挂载 Flags)
        VkDescriptorSetLayoutCreateInfo layoutInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO};
        layoutInfo.bindingCount = static_cast<uint32_t>(m_bindings.size());
        layoutInfo.pBindings = m_bindings.data();

        VkDescriptorSetLayoutBindingFlagsCreateInfo bindingFlagsInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_BINDING_FLAGS_CREATE_INFO};
        if (m_useBindless) {
            bindingFlagsInfo.bindingCount = static_cast<uint32_t>(m_bindingFlags.size());
            bindingFlagsInfo.pBindingFlags = m_bindingFlags.data();
            layoutInfo.pNext = &bindingFlagsInfo; // ★ 挂载 Bindless 标志
        }

        if (vkCreateDescriptorSetLayout(device->GetNativeDevice(), &layoutInfo, nullptr, &outLayout) != VK_SUCCESS) {
            throw std::runtime_error("Failed to create descriptor set layout via Builder!");
        }

        // 2. 分配 Set
        VkDescriptorSetAllocateInfo allocInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO};
        allocInfo.descriptorPool = pool;
        allocInfo.descriptorSetCount = 1;
        allocInfo.pSetLayouts = &outLayout;

        if (vkAllocateDescriptorSets(device->GetNativeDevice(), &allocInfo, &outSet) != VK_SUCCESS) {
            throw std::runtime_error("Failed to allocate descriptor set via Builder!");
        }

        // 3. 填装数据：由于分配完才知道 dstSet，所以在这里统一赋值并写入
        for (auto& write : m_writes) {
            write.dstSet = outSet;
        }

        vkUpdateDescriptorSets(device->GetNativeDevice(), static_cast<uint32_t>(m_writes.size()), m_writes.data(), 0, nullptr);

        return true;
    }

} // namespace Lizeral
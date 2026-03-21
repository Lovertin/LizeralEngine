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

        VkDescriptorSetLayoutBinding layoutBinding{};
        layoutBinding.binding = binding;
        layoutBinding.descriptorCount = layoutCount;
        layoutBinding.descriptorType = type;
        layoutBinding.stageFlags = stageFlags;
        m_bindings.push_back(layoutBinding);

        if (isBindless) {
            m_bindingFlags.push_back(VK_DESCRIPTOR_BINDING_PARTIALLY_BOUND_BIT);
            m_useBindless = true;
        } else {
            m_bindingFlags.push_back(0); 
        }
        VkWriteDescriptorSet write{VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET};
        write.dstBinding = binding;
        write.dstArrayElement = 0;
        write.descriptorType = type;
        write.descriptorCount = writeCount; 
        write.pImageInfo = imageInfos;
        m_writes.push_back(write);

        return *this;
    }

    VulkanDescriptorBuilder& VulkanDescriptorBuilder::BindImage(
        uint32_t binding, const VkDescriptorImageInfo* imageInfo, 
        VkDescriptorType type, VkShaderStageFlags stageFlags) 
    {
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
        std::unordered_map<VkDescriptorType, uint32_t> poolSizeMap;
        for (const auto& binding : m_bindings) {
            poolSizeMap[binding.descriptorType] += binding.descriptorCount;
        }

        std::vector<VkDescriptorPoolSize> poolSizes;
        for (const auto& pair : poolSizeMap) {
            poolSizes.push_back({ pair.first, pair.second });
        }

        VkDescriptorPoolCreateInfo poolInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO};
        poolInfo.poolSizeCount = static_cast<uint32_t>(poolSizes.size());
        poolInfo.pPoolSizes = poolSizes.data();
        poolInfo.maxSets = 1; 

        if (vkCreateDescriptorPool(device->GetNativeDevice(), &poolInfo, nullptr, &outPool) != VK_SUCCESS) {
            throw std::runtime_error("Failed to create descriptor pool via Builder!");
        }

        return Build(device, outPool, outLayout, outSet);
    }

    bool VulkanDescriptorBuilder::Build(VulkanDevice* device, VkDescriptorPool pool, VkDescriptorSetLayout& outLayout, VkDescriptorSet& outSet) {
        VkDescriptorSetLayoutCreateInfo layoutInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO};
        layoutInfo.bindingCount = static_cast<uint32_t>(m_bindings.size());
        layoutInfo.pBindings = m_bindings.data();

        VkDescriptorSetLayoutBindingFlagsCreateInfo bindingFlagsInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_BINDING_FLAGS_CREATE_INFO};
        if (m_useBindless) {
            bindingFlagsInfo.bindingCount = static_cast<uint32_t>(m_bindingFlags.size());
            bindingFlagsInfo.pBindingFlags = m_bindingFlags.data();
            layoutInfo.pNext = &bindingFlagsInfo;
        }

        if (vkCreateDescriptorSetLayout(device->GetNativeDevice(), &layoutInfo, nullptr, &outLayout) != VK_SUCCESS) {
            throw std::runtime_error("Failed to create descriptor set layout via Builder!");
        }

        VkDescriptorSetAllocateInfo allocInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO};
        allocInfo.descriptorPool = pool;
        allocInfo.descriptorSetCount = 1;
        allocInfo.pSetLayouts = &outLayout;

        if (vkAllocateDescriptorSets(device->GetNativeDevice(), &allocInfo, &outSet) != VK_SUCCESS) {
            throw std::runtime_error("Failed to allocate descriptor set via Builder!");
        }

        for (auto& write : m_writes) {
            write.dstSet = outSet;
        }

        vkUpdateDescriptorSets(device->GetNativeDevice(), static_cast<uint32_t>(m_writes.size()), m_writes.data(), 0, nullptr);

        return true;
    }

    VulkanDescriptorBuilder& VulkanDescriptorBuilder::BindAccelerationStructure(
        uint32_t binding, const VkAccelerationStructureKHR* as, VkShaderStageFlags stageFlags) 
    {
        VkDescriptorSetLayoutBinding layoutBinding{};
        layoutBinding.binding = binding;
        layoutBinding.descriptorCount = 1;
        layoutBinding.descriptorType = VK_DESCRIPTOR_TYPE_ACCELERATION_STRUCTURE_KHR;
        layoutBinding.stageFlags = stageFlags;
        m_bindings.push_back(layoutBinding);
        m_bindingFlags.push_back(0); 

        VkWriteDescriptorSetAccelerationStructureKHR asWrite{};
        asWrite.sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET_ACCELERATION_STRUCTURE_KHR;
        asWrite.accelerationStructureCount = 1;
        asWrite.pAccelerationStructures = as; 
        m_asWrites.push_back(asWrite);

        VkWriteDescriptorSet write{VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET};
        write.dstBinding = binding;
        write.dstArrayElement = 0;
        write.descriptorType = VK_DESCRIPTOR_TYPE_ACCELERATION_STRUCTURE_KHR;
        write.descriptorCount = 1;
        write.pNext = &m_asWrites.back();
        m_writes.push_back(write);

        return *this;
    }

} // namespace Lizeral
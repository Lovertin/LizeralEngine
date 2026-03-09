#include "VulkanComputePipeline.h"
#include "VulkanDevice.h"
#include <fstream>
#include <stdexcept>
#include <iostream>

namespace Lizeral {

    VulkanComputePipeline::VulkanComputePipeline(VulkanDevice* device, const std::string& shaderPath) 
        : m_device(device) {
        
        // 1. 读取并创建 Shader 模块
        auto shaderCode = readFile(shaderPath);
        VkShaderModule computeShaderModule = createShaderModule(shaderCode);

        VkPipelineShaderStageCreateInfo shaderStageInfo{};
        shaderStageInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
        shaderStageInfo.stage = VK_SHADER_STAGE_COMPUTE_BIT;
        shaderStageInfo.module = computeShaderModule;
        shaderStageInfo.pName = "main"; // 对应 glsl 里的 void main()

        // 2. 核心魔法：定义 Push Constant，用来接收 64位的 BDA 地址
        VkPushConstantRange pushConstantRange{};
        pushConstantRange.stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;
        pushConstantRange.offset = 0;
        pushConstantRange.size = sizeof(uint64_t); // 只需要 8 个字节！

        // 3. 创建 Pipeline Layout (管线布局)
        VkPipelineLayoutCreateInfo pipelineLayoutInfo{};
        pipelineLayoutInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO;
        pipelineLayoutInfo.pushConstantRangeCount = 1;
        pipelineLayoutInfo.pPushConstantRanges = &pushConstantRange;
        // 注意：这里 pSetLayouts 是空的！没有任何 Descriptor 绑定！
        pipelineLayoutInfo.setLayoutCount = 0;
        pipelineLayoutInfo.pSetLayouts = nullptr;

        if (vkCreatePipelineLayout(m_device->GetNativeDevice(), &pipelineLayoutInfo, nullptr, &m_pipelineLayout) != VK_SUCCESS) {
            throw std::runtime_error("Failed to create compute pipeline layout!");
        }

        // 4. 创建最终的 Compute Pipeline
        VkComputePipelineCreateInfo pipelineInfo{};
        pipelineInfo.sType = VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO;
        pipelineInfo.layout = m_pipelineLayout;
        pipelineInfo.stage = shaderStageInfo;

        if (vkCreateComputePipelines(m_device->GetNativeDevice(), VK_NULL_HANDLE, 1, &pipelineInfo, nullptr, &m_pipeline) != VK_SUCCESS) {
            throw std::runtime_error("Failed to create compute pipeline!");
        }

        // Shader 模块在管线创建完后就可以销毁了
        vkDestroyShaderModule(m_device->GetNativeDevice(), computeShaderModule, nullptr);
        std::cout << "[VulkanCompute] BDA Compute Pipeline created successfully!" << std::endl;
    }

    VulkanComputePipeline::~VulkanComputePipeline() {
        if (m_pipeline != VK_NULL_HANDLE) {
            vkDestroyPipeline(m_device->GetNativeDevice(), m_pipeline, nullptr);
        }
        if (m_pipelineLayout != VK_NULL_HANDLE) {
            vkDestroyPipelineLayout(m_device->GetNativeDevice(), m_pipelineLayout, nullptr);
        }
    }

    std::vector<char> VulkanComputePipeline::readFile(const std::string& filename) {
        std::ifstream file(filename, std::ios::ate | std::ios::binary);
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open shader file: " + filename);
        }
        size_t fileSize = (size_t)file.tellg();
        std::vector<char> buffer(fileSize);
        file.seekg(0);
        file.read(buffer.data(), fileSize);
        file.close();
        return buffer;
    }

    VkShaderModule VulkanComputePipeline::createShaderModule(const std::vector<char>& code) {
        VkShaderModuleCreateInfo createInfo{};
        createInfo.sType = VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO;
        createInfo.codeSize = code.size();
        createInfo.pCode = reinterpret_cast<const uint32_t*>(code.data());

        VkShaderModule shaderModule;
        if (vkCreateShaderModule(m_device->GetNativeDevice(), &createInfo, nullptr, &shaderModule) != VK_SUCCESS) {
            throw std::runtime_error("Failed to create shader module!");
        }
        return shaderModule;
    }

} // namespace Lizeral
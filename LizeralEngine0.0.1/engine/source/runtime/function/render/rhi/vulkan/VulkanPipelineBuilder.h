#pragma once
#include <vulkan/vulkan.h>
#include <vector>

namespace Lizeral {

    class VulkanDevice;

    class VulkanPipelineBuilder {
    public:
        // 构造时，我们会给所有状态填上“现代引擎最常用的默认值”
        VulkanPipelineBuilder();

        // --- 核心配置接口 (链式调用) ---

        // 1. 添加 Shader (可以加 Vert/Frag，以后也能加 Task/Mesh)
        VulkanPipelineBuilder& AddShaderStage(VkShaderStageFlagBits stage, VkShaderModule shaderModule);

        // 2. 传统顶点输入 (对于我们的无Buffer测试和未来的Mesh Shader，我们根本不调它！)
        VulkanPipelineBuilder& SetVertexInput(const VkPipelineVertexInputStateCreateInfo& vertexInputInfo);

        // 3. 图元装配 (点、线、三角形)
        VulkanPipelineBuilder& SetInputAssembly(VkPrimitiveTopology topology, bool primitiveRestartEnable = false);

        // 4. 光栅化设置 (线框/填充，背面剔除)
        VulkanPipelineBuilder& SetRasterization(VkPolygonMode polygonMode, VkCullModeFlags cullMode, VkFrontFace frontFace = VK_FRONT_FACE_COUNTER_CLOCKWISE);

        // 5. 多重采样 (抗锯齿，默认关)
        VulkanPipelineBuilder& SetMultisampling(VkSampleCountFlagBits samples = VK_SAMPLE_COUNT_1_BIT);

        // 6. 深度与模板测试
        VulkanPipelineBuilder& SetDepthStencil(bool depthTestEnable, bool depthWriteEnable, VkCompareOp depthCompareOp);

        // 7. 颜色混合 (透明度处理)
        VulkanPipelineBuilder& SetColorBlend(bool blendEnable, VkBlendOp colorBlendOp = VK_BLEND_OP_ADD, 
                                             VkBlendFactor srcColorBlendFactor = VK_BLEND_FACTOR_SRC_ALPHA, 
                                             VkBlendFactor dstColorBlendFactor = VK_BLEND_FACTOR_ONE_MINUS_SRC_ALPHA);

        // 8. 绑定管线布局 (非常重要：用来传 Push Constants 和 Descriptor Sets)
        VulkanPipelineBuilder& SetPipelineLayout(VkPipelineLayout pipelineLayout);

        // --- 终极兵工厂产出 ---
        // 传入 Device 和 RenderPass，正式将所有状态“焊死”成一条物理硬件管线！
        VkPipeline Build(VulkanDevice* device, VkFormat colorFormat, VkFormat depthFormat); //dynamic type

    private:
        std::vector<VkPipelineShaderStageCreateInfo> m_shaderStages;
        
        // 各种状态机的结构体缓存
        VkPipelineVertexInputStateCreateInfo m_vertexInputInfo;
        VkPipelineInputAssemblyStateCreateInfo m_inputAssembly;
        VkPipelineViewportStateCreateInfo m_viewportState;
        VkPipelineRasterizationStateCreateInfo m_rasterizer;
        VkPipelineMultisampleStateCreateInfo m_multisampling;
        VkPipelineDepthStencilStateCreateInfo m_depthStencil;
        VkPipelineColorBlendAttachmentState m_colorBlendAttachment;
        VkPipelineColorBlendStateCreateInfo m_colorBlending;
        
        // 动态状态 (现代引擎标配)
        std::vector<VkDynamicState> m_dynamicStates;
        VkPipelineDynamicStateCreateInfo m_dynamicStateInfo;

        VkPipelineLayout m_pipelineLayout { VK_NULL_HANDLE };
    };

} // namespace Lizeral
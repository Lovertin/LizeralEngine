#include "VulkanRenderPass.h"
#include "VulkanDevice.h"
#include "VulkanSwapchain.h"
#include <stdexcept>
#include <iostream>

namespace Lizeral {

    VulkanRenderPass::VulkanRenderPass(VulkanDevice* device, VulkanSwapchain* swapchain) 
        : m_device(device) {
        
        // 1. 定义颜色附件 (我们要画的靶子)
        VkAttachmentDescription colorAttachment{};
        colorAttachment.format = swapchain->GetImageFormat(); // 必须和 Swapchain 格式一致
        colorAttachment.samples = VK_SAMPLE_COUNT_1_BIT;      // 不开多重采样抗锯齿 (MSAA)
        
        // 核心理论落地：
        colorAttachment.loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR;   // 绘制开始前，把这块显存清空（通常清空为黑色）
        colorAttachment.storeOp = VK_ATTACHMENT_STORE_OP_STORE; // 绘制结束后，把结果保存在显存里，为了给显示器看
        
        // 模板测试(Stencil)相关的操作，我们现在不用，选 DONT_CARE 提升性能
        colorAttachment.stencilLoadOp = VK_ATTACHMENT_LOAD_OP_DONT_CARE;
        colorAttachment.stencilStoreOp = VK_ATTACHMENT_STORE_OP_DONT_CARE;
        
        // 内存布局转换的魔法：
        colorAttachment.initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;      // 刚拿过来时，我们不在乎里面原本是什么乱码
        colorAttachment.finalLayout = VK_IMAGE_LAYOUT_PRESENT_SRC_KHR;  // 画完后，自动变形成“准备送给屏幕显示的布局”

        // 2. 声明附件的引用 (在 Subpass 中使用)
        VkAttachmentReference colorAttachmentRef{};
        colorAttachmentRef.attachment = 0; // 对应上面 colorAttachment 在数组里的索引 (目前只有 1 个)
        colorAttachmentRef.layout = VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL; // 在画图期间，要求它处于“最佳涂色布局”

        // 3. 定义子通道 (Subpass)
        VkSubpassDescription subpass{};
        subpass.pipelineBindPoint = VK_PIPELINE_BIND_POINT_GRAPHICS; // 这是一个图形渲染通道，不是 Compute
        subpass.colorAttachmentCount = 1;
        subpass.pColorAttachments = &colorAttachmentRef; // 把靶子挂上去

        // 4. 定义子通道依赖 (Subpass Dependency) - 极其关键的同步机制！
        // 显卡在后台准备图片时需要时间，我们不能在图片还没准备好的时候就往上画。
        VkSubpassDependency dependency{};
        dependency.srcSubpass = VK_SUBPASS_EXTERNAL; // 代表 RenderPass 开始之前的那个“外部隐形操作”
        dependency.dstSubpass = 0;                   // 代表我们自己的这个 Subpass
        
        // 我们要等“读取 Swapchain 图像”这个操作完成
        dependency.srcStageMask = VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT;
        dependency.srcAccessMask = 0;
        // 然后我们才能进行“写颜色”操作
        dependency.dstStageMask = VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT;
        dependency.dstAccessMask = VK_ACCESS_COLOR_ATTACHMENT_WRITE_BIT;

        // 5. 拼装出完整的 Render Pass！
        VkRenderPassCreateInfo renderPassInfo{};
        renderPassInfo.sType = VK_STRUCTURE_TYPE_RENDER_PASS_CREATE_INFO;
        renderPassInfo.attachmentCount = 1;
        renderPassInfo.pAttachments = &colorAttachment;
        renderPassInfo.subpassCount = 1;
        renderPassInfo.pSubpasses = &subpass;
        renderPassInfo.dependencyCount = 1;
        renderPassInfo.pDependencies = &dependency;

        if (vkCreateRenderPass(m_device->GetNativeDevice(), &renderPassInfo, nullptr, &m_renderPass) != VK_SUCCESS) {
            throw std::runtime_error("Failed to create render pass!");
        }

        std::cout << "[VulkanRenderPass] Default Swapchain Render Pass created successfully!" << std::endl;
    }

    VulkanRenderPass::~VulkanRenderPass() {
        if (m_renderPass != VK_NULL_HANDLE) {
            vkDestroyRenderPass(m_device->GetNativeDevice(), m_renderPass, nullptr);
        }
    }

} // namespace Lizeral
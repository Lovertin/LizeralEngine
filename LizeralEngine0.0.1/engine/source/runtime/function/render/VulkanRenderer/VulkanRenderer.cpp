#include "VulkanRenderer.h"
#include "runtime/function/render/rhi/vulkan/VulkanContext.h"
#include "runtime/function/render/rhi/vulkan/VulkanDevice.h"
#include "runtime/function/render/rhi/vulkan/VulkanCommandPool.h"
#include "runtime/function/render/rhi/vulkan/VulkanCommandBuffer.h"
#include "runtime/function/render/rhi/vulkan/VulkanSwapchain.h"

#include <stdexcept>
#include <iostream>
#include <array>

namespace Lizeral {

    VulkanRenderer::VulkanRenderer(VulkanContext* context, VulkanDevice* device, GLFWwindow* window)
        : m_context(context), m_device(device), m_window(window) {
        
        RecreateSwapchain();
        CreateCommandBuffers();
        CreateSyncObjects();

        std::cout << "[VulkanRenderer] Dynamic Renderer Initialized Successfully!" << std::endl;
    }

    VulkanRenderer::~VulkanRenderer() {
        vkDeviceWaitIdle(m_device->GetNativeDevice());

        FreeCommandBuffers();
        CleanupDepthResources();

        for (int i = 0; i < MAX_FRAMES_IN_FLIGHT; i++) {
            vkDestroySemaphore(m_device->GetNativeDevice(), m_imageAvailableSemaphores[i], nullptr);
            vkDestroyFence(m_device->GetNativeDevice(), m_inFlightFences[i], nullptr);
        }
        for (auto& semaphore : m_renderFinishedSemaphores) {
            vkDestroySemaphore(m_device->GetNativeDevice(), semaphore, nullptr);
        }
    }

    VkFormat VulkanRenderer::GetSwapchainFormat() const {
        return m_swapchain->GetImageFormat();
    }

    VkExtent2D VulkanRenderer::GetSwapchainExtent() const {
        return m_swapchain->GetExtent();
    }

    void VulkanRenderer::RecreateSwapchain(int width, int height) {
        if (width == 0 || height == 0) return; 

        vkDeviceWaitIdle(m_device->GetNativeDevice());

        m_swapchain.reset();
        CleanupDepthResources();

        m_swapchain = std::make_unique<VulkanSwapchain>(m_context, m_device, m_device->GetSurface(), width, height);

        // ★ 关键修复：向 Swapchain 索要它最终决定的真实物理分辨率
        VkExtent2D actualExtent = m_swapchain->GetExtent();

        VkImageCreateInfo imageInfo{};
        imageInfo.sType = VK_STRUCTURE_TYPE_IMAGE_CREATE_INFO;
        imageInfo.imageType = VK_IMAGE_TYPE_2D;
        // ★ 使用真实的 actualExtent，绝不能用外部传进来的 width/height
        imageInfo.extent.width = actualExtent.width;
        imageInfo.extent.height = actualExtent.height;
        imageInfo.extent.depth = 1;
        imageInfo.mipLevels = 1;
        imageInfo.arrayLayers = 1;
        imageInfo.format = VK_FORMAT_D32_SFLOAT; 
        imageInfo.tiling = VK_IMAGE_TILING_OPTIMAL;
        imageInfo.initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;
        imageInfo.usage = VK_IMAGE_USAGE_DEPTH_STENCIL_ATTACHMENT_BIT;
        imageInfo.samples = VK_SAMPLE_COUNT_1_BIT;
        imageInfo.sharingMode = VK_SHARING_MODE_EXCLUSIVE;
        vkCreateImage(m_device->GetNativeDevice(), &imageInfo, nullptr, &m_depthImage);

        VkMemoryRequirements memRequirements;
        vkGetImageMemoryRequirements(m_device->GetNativeDevice(), m_depthImage, &memRequirements);
        
        VkPhysicalDeviceMemoryProperties memProperties;
        vkGetPhysicalDeviceMemoryProperties(m_context->GetPhysicalDevice(), &memProperties);
        uint32_t memTypeIndex = 0;
        for (uint32_t i = 0; i < memProperties.memoryTypeCount; i++) {
            if ((memRequirements.memoryTypeBits & (1 << i)) && 
                (memProperties.memoryTypes[i].propertyFlags & VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT) == VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT) {
                memTypeIndex = i; break;
            }
        }
        VkMemoryAllocateInfo allocInfo{VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO};
        allocInfo.allocationSize = memRequirements.size;
        allocInfo.memoryTypeIndex = memTypeIndex;
        vkAllocateMemory(m_device->GetNativeDevice(), &allocInfo, nullptr, &m_depthImageMemory);
        vkBindImageMemory(m_device->GetNativeDevice(), m_depthImage, m_depthImageMemory, 0);

        VkImageViewCreateInfo viewInfo{};
        viewInfo.sType = VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO;
        viewInfo.image = m_depthImage;
        viewInfo.viewType = VK_IMAGE_VIEW_TYPE_2D;
        viewInfo.format = VK_FORMAT_D32_SFLOAT;
        viewInfo.subresourceRange.aspectMask = VK_IMAGE_ASPECT_DEPTH_BIT;
        viewInfo.subresourceRange.baseMipLevel = 0;
        viewInfo.subresourceRange.levelCount = 1;
        viewInfo.subresourceRange.baseArrayLayer = 0;
        viewInfo.subresourceRange.layerCount = 1;
        vkCreateImageView(m_device->GetNativeDevice(), &viewInfo, nullptr, &m_depthImageView);
        
        // ★ 彻底结束！再也没有 Framebuffer 和 RenderPass 了！
    }

    void VulkanRenderer::CleanupDepthResources() {
        if (m_depthImageView != VK_NULL_HANDLE) {
            vkDestroyImageView(m_device->GetNativeDevice(), m_depthImageView, nullptr);
            vkDestroyImage(m_device->GetNativeDevice(), m_depthImage, nullptr);
            vkFreeMemory(m_device->GetNativeDevice(), m_depthImageMemory, nullptr);
            m_depthImageView = VK_NULL_HANDLE;
            m_depthImage = VK_NULL_HANDLE;
            m_depthImageMemory = VK_NULL_HANDLE;
        }
    }

    void VulkanRenderer::CreateCommandBuffers() {
        m_commandPool.reset(new VulkanCommandPool(m_device));
        m_commandBuffers.resize(MAX_FRAMES_IN_FLIGHT);
        for (int i = 0; i < MAX_FRAMES_IN_FLIGHT; i++) {
            m_commandBuffers[i] = std::make_unique<VulkanCommandBuffer>(m_device, m_commandPool.get());
        }
    }

    void VulkanRenderer::FreeCommandBuffers() {
        m_commandBuffers.clear(); 
        m_commandPool.reset();
    }

    void VulkanRenderer::CreateSyncObjects() {
        m_imageAvailableSemaphores.resize(MAX_FRAMES_IN_FLIGHT);
        m_inFlightFences.resize(MAX_FRAMES_IN_FLIGHT);
        // 初始化时，如果还没创建 Swapchain，先随便给个值，真正渲染前会保证其正确
        size_t imageCount = m_swapchain ? m_swapchain->GetImageViews().size() : MAX_FRAMES_IN_FLIGHT;
        m_renderFinishedSemaphores.resize(imageCount);

        VkSemaphoreCreateInfo semaphoreInfo{};
        semaphoreInfo.sType = VK_STRUCTURE_TYPE_SEMAPHORE_CREATE_INFO;

        VkFenceCreateInfo fenceInfo{};
        fenceInfo.sType = VK_STRUCTURE_TYPE_FENCE_CREATE_INFO;
        fenceInfo.flags = VK_FENCE_CREATE_SIGNALED_BIT; 

        for (int i = 0; i < MAX_FRAMES_IN_FLIGHT; i++) {
            if (vkCreateSemaphore(m_device->GetNativeDevice(), &semaphoreInfo, nullptr, &m_imageAvailableSemaphores[i]) != VK_SUCCESS ||
                vkCreateFence(m_device->GetNativeDevice(), &fenceInfo, nullptr, &m_inFlightFences[i]) != VK_SUCCESS) {
                throw std::runtime_error("Failed to create synchronization objects for a frame!");
            }
        }
        for (size_t i = 0; i < m_renderFinishedSemaphores.size(); i++) {
            if (vkCreateSemaphore(m_device->GetNativeDevice(), &semaphoreInfo, nullptr, &m_renderFinishedSemaphores[i]) != VK_SUCCESS) {
                throw std::runtime_error("Failed to create render finished semaphores!");
            }
        }
    }

    VkCommandBuffer VulkanRenderer::GetCurrentCommandBuffer() const {
        if (!m_isFrameStarted) throw std::runtime_error("Cannot get command buffer when frame not in progress!");
        return m_commandBuffers[m_currentFrame]->GetNativeBuffer();
    }

    VkCommandBuffer VulkanRenderer::BeginFrame() {
        if (m_isFrameStarted) throw std::runtime_error("Can't call BeginFrame while frame is already in progress!");

        vkWaitForFences(m_device->GetNativeDevice(), 1, &m_inFlightFences[m_currentFrame], VK_TRUE, UINT64_MAX);

        VkResult result = vkAcquireNextImageKHR(m_device->GetNativeDevice(), m_swapchain->GetNativeSwapchain(), UINT64_MAX, m_imageAvailableSemaphores[m_currentFrame], VK_NULL_HANDLE, &m_imageIndex);

        if (result == VK_ERROR_OUT_OF_DATE_KHR) {
            RecreateSwapchain();
            return VK_NULL_HANDLE; 
        } else if (result != VK_SUCCESS && result != VK_SUBOPTIMAL_KHR) {
            throw std::runtime_error("Failed to acquire swapchain image!");
        }

        vkResetFences(m_device->GetNativeDevice(), 1, &m_inFlightFences[m_currentFrame]);

        m_isFrameStarted = true;
        auto commandBuffer = GetCurrentCommandBuffer();
        vkResetCommandBuffer(commandBuffer, 0);

        VkCommandBufferBeginInfo beginInfo{};
        beginInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
        if (vkBeginCommandBuffer(commandBuffer, &beginInfo) != VK_SUCCESS) {
            throw std::runtime_error("Failed to begin recording command buffer!");
        }

        return commandBuffer;
    }

    void VulkanRenderer::EndFrame() {
        if (!m_isFrameStarted) throw std::runtime_error("Can't call EndFrame while frame is not in progress!");

        auto commandBuffer = GetCurrentCommandBuffer();
        if (vkEndCommandBuffer(commandBuffer) != VK_SUCCESS) {
            throw std::runtime_error("Failed to record command buffer!");
        }

        VkSubmitInfo submitInfo{};
        submitInfo.sType = VK_STRUCTURE_TYPE_SUBMIT_INFO;

        VkSemaphore waitSemaphores[] = { m_imageAvailableSemaphores[m_currentFrame] };
        VkPipelineStageFlags waitStages[] = { VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT };
        submitInfo.waitSemaphoreCount = 1;
        submitInfo.pWaitSemaphores = waitSemaphores;
        submitInfo.pWaitDstStageMask = waitStages;

        submitInfo.commandBufferCount = 1;
        submitInfo.pCommandBuffers = &commandBuffer;

        VkSemaphore signalSemaphores[] = { m_renderFinishedSemaphores[m_imageIndex] };
        submitInfo.signalSemaphoreCount = 1;
        submitInfo.pSignalSemaphores = signalSemaphores;

        if (vkQueueSubmit(m_device->GetGraphicsQueue(), 1, &submitInfo, m_inFlightFences[m_currentFrame]) != VK_SUCCESS) {
            throw std::runtime_error("Failed to submit draw command buffer!");
        }

        VkPresentInfoKHR presentInfo{};
        presentInfo.sType = VK_STRUCTURE_TYPE_PRESENT_INFO_KHR;
        presentInfo.waitSemaphoreCount = 1;
        presentInfo.pWaitSemaphores = signalSemaphores;

        VkSwapchainKHR swapchains[] = { m_swapchain->GetNativeSwapchain() };
        presentInfo.swapchainCount = 1;
        presentInfo.pSwapchains = swapchains;
        presentInfo.pImageIndices = &m_imageIndex;

        VkResult result = vkQueuePresentKHR(m_device->GetPresentQueue(), &presentInfo);

        if (result == VK_ERROR_OUT_OF_DATE_KHR || result == VK_SUBOPTIMAL_KHR) {
            RecreateSwapchain();
        } else if (result != VK_SUCCESS) {
            throw std::runtime_error("Failed to present swapchain image!");
        }

        m_isFrameStarted = false;
        m_currentFrame = (m_currentFrame + 1) % MAX_FRAMES_IN_FLIGHT;
    }

    // =========================================================================
    // ★ 动态渲染核心：管线状态转移与绘制开始
    // =========================================================================
    void VulkanRenderer::BeginRendering(VkCommandBuffer commandBuffer) {
        if (!m_isFrameStarted || commandBuffer != GetCurrentCommandBuffer()) {
            throw std::runtime_error("Can't begin rendering on an invalid command buffer!");
        }

        // 1. 手动转换内存布局 (Layout Transition) - 使用 Vulkan 1.3 同步 2.0
        VkImageMemoryBarrier2 colorBarrier{};
        colorBarrier.sType = VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER_2;
        colorBarrier.srcStageMask = VK_PIPELINE_STAGE_2_COLOR_ATTACHMENT_OUTPUT_BIT;
        colorBarrier.srcAccessMask = 0;
        colorBarrier.dstStageMask = VK_PIPELINE_STAGE_2_COLOR_ATTACHMENT_OUTPUT_BIT;
        colorBarrier.dstAccessMask = VK_ACCESS_2_COLOR_ATTACHMENT_WRITE_BIT;
        colorBarrier.oldLayout = m_swapchain->GetImageLayout(m_imageIndex); 
        colorBarrier.newLayout = VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL;
        colorBarrier.image = m_swapchain->GetNativeImages()[m_imageIndex]; 
        colorBarrier.subresourceRange = {VK_IMAGE_ASPECT_COLOR_BIT, 0, 1, 0, 1};

        VkImageMemoryBarrier2 depthBarrier{};
        depthBarrier.sType = VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER_2;
        depthBarrier.srcStageMask = VK_PIPELINE_STAGE_2_EARLY_FRAGMENT_TESTS_BIT | VK_PIPELINE_STAGE_2_LATE_FRAGMENT_TESTS_BIT;
        depthBarrier.srcAccessMask = 0;
        depthBarrier.dstStageMask = VK_PIPELINE_STAGE_2_EARLY_FRAGMENT_TESTS_BIT | VK_PIPELINE_STAGE_2_LATE_FRAGMENT_TESTS_BIT;
        depthBarrier.dstAccessMask = VK_ACCESS_2_DEPTH_STENCIL_ATTACHMENT_WRITE_BIT;
        depthBarrier.oldLayout = VK_IMAGE_LAYOUT_UNDEFINED;
        depthBarrier.newLayout = VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL;
        depthBarrier.image = m_depthImage;
        depthBarrier.subresourceRange = {VK_IMAGE_ASPECT_DEPTH_BIT, 0, 1, 0, 1};

        VkImageMemoryBarrier2 barriers[] = {colorBarrier, depthBarrier};
        VkDependencyInfo depInfo{};
        depInfo.sType = VK_STRUCTURE_TYPE_DEPENDENCY_INFO;
        depInfo.imageMemoryBarrierCount = 2;
        depInfo.pImageMemoryBarriers = barriers;
        vkCmdPipelineBarrier2(commandBuffer, &depInfo); // 瞬间完成变形！

        // 2. 组装临时画板并开画！
        VkRenderingAttachmentInfo colorAttachment{};
        colorAttachment.sType = VK_STRUCTURE_TYPE_RENDERING_ATTACHMENT_INFO;
        colorAttachment.imageView = m_swapchain->GetImageViews()[m_imageIndex];
        colorAttachment.imageLayout = VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL;
        colorAttachment.loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR;
        colorAttachment.storeOp = VK_ATTACHMENT_STORE_OP_STORE;
        colorAttachment.clearValue.color = {{0.05f, 0.05f, 0.05f, 1.0f}};

        VkRenderingAttachmentInfo depthAttachment{};
        depthAttachment.sType = VK_STRUCTURE_TYPE_RENDERING_ATTACHMENT_INFO;
        depthAttachment.imageView = m_depthImageView;
        depthAttachment.imageLayout = VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL;
        depthAttachment.loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR;
        depthAttachment.storeOp = VK_ATTACHMENT_STORE_OP_DONT_CARE;
        depthAttachment.clearValue.depthStencil = {0.0f, 0};

        VkRenderingInfo renderInfo{};
        renderInfo.sType = VK_STRUCTURE_TYPE_RENDERING_INFO;
        renderInfo.renderArea.offset = {0, 0};
        renderInfo.renderArea.extent = m_swapchain->GetExtent();
        renderInfo.layerCount = 1;
        renderInfo.colorAttachmentCount = 1;
        renderInfo.pColorAttachments = &colorAttachment;
        renderInfo.pDepthAttachment = &depthAttachment;

        vkCmdBeginRendering(commandBuffer, &renderInfo);

        // 设置视口
        VkViewport viewport{0.0f, 0.0f, static_cast<float>(m_swapchain->GetExtent().width), static_cast<float>(m_swapchain->GetExtent().height), 0.0f, 1.0f};
        VkRect2D scissor{{0, 0}, m_swapchain->GetExtent()};
        vkCmdSetViewport(commandBuffer, 0, 1, &viewport);
        vkCmdSetScissor(commandBuffer, 0, 1, &scissor);
    }

    void VulkanRenderer::EndRendering(VkCommandBuffer commandBuffer) {
        if (!m_isFrameStarted || commandBuffer != GetCurrentCommandBuffer()) {
            throw std::runtime_error("Can't end rendering on an invalid command buffer!");
        }

        vkCmdEndRendering(commandBuffer);

        // 画完了，手动把颜色图片变形为“展示模式 (PRESENT)”
        VkImageMemoryBarrier2 presentBarrier{};
        presentBarrier.sType = VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER_2;
        presentBarrier.srcStageMask = VK_PIPELINE_STAGE_2_COLOR_ATTACHMENT_OUTPUT_BIT;
        presentBarrier.srcAccessMask = VK_ACCESS_2_COLOR_ATTACHMENT_WRITE_BIT;
        presentBarrier.dstStageMask = VK_PIPELINE_STAGE_2_BOTTOM_OF_PIPE_BIT;
        presentBarrier.dstAccessMask = 0;
        presentBarrier.oldLayout = VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL;
        presentBarrier.newLayout = VK_IMAGE_LAYOUT_PRESENT_SRC_KHR;
        presentBarrier.image = m_swapchain->GetNativeImages()[m_imageIndex]; 
        presentBarrier.subresourceRange = {VK_IMAGE_ASPECT_COLOR_BIT, 0, 1, 0, 1};

        VkDependencyInfo depInfo{};
        depInfo.sType = VK_STRUCTURE_TYPE_DEPENDENCY_INFO;
        depInfo.imageMemoryBarrierCount = 1;
        depInfo.pImageMemoryBarriers = &presentBarrier;
        vkCmdPipelineBarrier2(commandBuffer, &depInfo);
        m_swapchain->SetImageLayout(m_imageIndex, VK_IMAGE_LAYOUT_PRESENT_SRC_KHR);
    }

} // namespace Lizeral
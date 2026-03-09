#include "VulkanRenderer.h"
#include "runtime/function/render/rhi/vulkan/VulkanContext.h"
#include "runtime/function/render/rhi/vulkan/VulkanDevice.h"
#include "runtime/function/render/rhi/vulkan/VulkanSwapchain.h"
#include "runtime/function/render/rhi/vulkan/VulkanRenderPass.h"
#include "runtime/function/render/rhi/vulkan/VulkanFramebuffer.h"
#include "runtime/function/render/rhi/vulkan/VulkanCommandPool.h"
#include "runtime/function/render/rhi/vulkan/VulkanCommandBuffer.h"

#include <stdexcept>
#include <iostream>
#include <array>

namespace Lizeral {

    VulkanRenderer::VulkanRenderer(VulkanContext* context, VulkanDevice* device, GLFWwindow* window)
        : m_context(context), m_device(device), m_window(window) {
        
        RecreateSwapchain();
        CreateCommandBuffers();
        CreateSyncObjects();

        std::cout << "[VulkanRenderer] Renderer Initialized Successfully!" << std::endl;
    }

    VulkanRenderer::~VulkanRenderer() {
        // 在销毁任何东西前，必须等 GPU 彻底停下手中的活
        vkDeviceWaitIdle(m_device->GetNativeDevice());

        FreeCommandBuffers();

        // 销毁同步对象
        for (int i = 0; i < MAX_FRAMES_IN_FLIGHT; i++) {
            vkDestroySemaphore(m_device->GetNativeDevice(), m_imageAvailableSemaphores[i], nullptr);
            vkDestroyFence(m_device->GetNativeDevice(), m_inFlightFences[i], nullptr);
        }
        for (auto& semaphore : m_renderFinishedSemaphores) {
            vkDestroySemaphore(m_device->GetNativeDevice(), semaphore, nullptr);
        }

        // unique_ptr 会在这一步自动按逆序安全释放 (Swapchain, RenderPass, Framebuffer, CommandPool)
    }

    void VulkanRenderer::RecreateSwapchain() {
        // 处理窗口最小化暂停渲染的逻辑
        int width = 0, height = 0;
        glfwGetFramebufferSize(m_window, &width, &height);
        while (width == 0 || height == 0) {
            glfwGetFramebufferSize(m_window, &width, &height);
            glfwWaitEvents();
        }

        vkDeviceWaitIdle(m_device->GetNativeDevice());

        m_defaultFramebuffers.reset();
        m_defaultRenderPass.reset();
        m_swapchain.reset();

        // 使用 unique_ptr 的 reset 来销毁旧对象并创建新对象
        m_swapchain.reset(new VulkanSwapchain(m_context, m_device, m_device->GetSurface(), width, height));
        m_defaultRenderPass.reset(new VulkanRenderPass(m_device, m_swapchain.get()));
        m_defaultFramebuffers.reset(new VulkanFramebuffer(m_device, m_defaultRenderPass.get(), m_swapchain.get()));
    }

    void VulkanRenderer::CreateCommandBuffers() {
        m_commandPool.reset(new VulkanCommandPool(m_device));
        m_commandBuffers.resize(MAX_FRAMES_IN_FLIGHT);
        for (int i = 0; i < MAX_FRAMES_IN_FLIGHT; i++) {
            m_commandBuffers[i] = std::make_unique<VulkanCommandBuffer>(m_device, m_commandPool.get());
        }
    }

    void VulkanRenderer::FreeCommandBuffers() {
        m_commandBuffers.clear(); // 触发 unique_ptr 析构
        m_commandPool.reset();
    }

    void VulkanRenderer::CreateSyncObjects() {
        m_imageAvailableSemaphores.resize(MAX_FRAMES_IN_FLIGHT);
        m_inFlightFences.resize(MAX_FRAMES_IN_FLIGHT);
        // 注意：呈现信号灯的数量与 Swapchain 里的图片数量严格一致！
        m_renderFinishedSemaphores.resize(m_swapchain->GetImageViews().size());

        VkSemaphoreCreateInfo semaphoreInfo{};
        semaphoreInfo.sType = VK_STRUCTURE_TYPE_SEMAPHORE_CREATE_INFO;

        VkFenceCreateInfo fenceInfo{};
        fenceInfo.sType = VK_STRUCTURE_TYPE_FENCE_CREATE_INFO;
        fenceInfo.flags = VK_FENCE_CREATE_SIGNALED_BIT; // 初始为绿灯

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

    VkCommandBuffer VulkanRenderer::BeginFrame() {
        if (m_isFrameStarted) {
            throw std::runtime_error("Can't call BeginFrame while frame is already in progress!");
        }

        // 1. 阻塞等待 GPU 干完上一帧的活
        vkWaitForFences(m_device->GetNativeDevice(), 1, &m_inFlightFences[m_currentFrame], VK_TRUE, UINT64_MAX);

        // 2. 问显示器要下一张画布
        VkResult result = vkAcquireNextImageKHR(m_device->GetNativeDevice(), m_swapchain->GetNativeSwapchain(), UINT64_MAX, m_imageAvailableSemaphores[m_currentFrame], VK_NULL_HANDLE, &m_imageIndex);

        // 处理窗口拖拽改变大小的情况
        if (result == VK_ERROR_OUT_OF_DATE_KHR) {
            RecreateSwapchain();
            return VK_NULL_HANDLE; // 告诉上层这一帧作废，重新开始
        } else if (result != VK_SUCCESS && result != VK_SUBOPTIMAL_KHR) {
            throw std::runtime_error("Failed to acquire swapchain image!");
        }

        // 3. 拿到图片了，关上铁门
        vkResetFences(m_device->GetNativeDevice(), 1, &m_inFlightFences[m_currentFrame]);

        // 4. 重置并开始录制 CommandBuffer
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
        if (!m_isFrameStarted) {
            throw std::runtime_error("Can't call EndFrame while frame is not in progress!");
        }

        auto commandBuffer = GetCurrentCommandBuffer();
        if (vkEndCommandBuffer(commandBuffer) != VK_SUCCESS) {
            throw std::runtime_error("Failed to record command buffer!");
        }

        // ==========================================================
        // ★ 核心提交流程
        // ==========================================================
        VkSubmitInfo submitInfo{};
        submitInfo.sType = VK_STRUCTURE_TYPE_SUBMIT_INFO;

        // 等待：图片准备好 -> 才能执行 颜色输出 阶段
        VkSemaphore waitSemaphores[] = { m_imageAvailableSemaphores[m_currentFrame] };
        VkPipelineStageFlags waitStages[] = { VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT };
        submitInfo.waitSemaphoreCount = 1;
        submitInfo.pWaitSemaphores = waitSemaphores;
        submitInfo.pWaitDstStageMask = waitStages;

        submitInfo.commandBufferCount = 1;
        submitInfo.pCommandBuffers = &commandBuffer;

        // 信号：画完了 -> 点亮 对应图片索引 的信号灯
        VkSemaphore signalSemaphores[] = { m_renderFinishedSemaphores[m_imageIndex] };
        submitInfo.signalSemaphoreCount = 1;
        submitInfo.pSignalSemaphores = signalSemaphores;

        if (vkQueueSubmit(m_device->GetGraphicsQueue(), 1, &submitInfo, m_inFlightFences[m_currentFrame]) != VK_SUCCESS) {
            throw std::runtime_error("Failed to submit draw command buffer!");
        }

        // ==========================================================
        // ★ 呈现给显示器
        // ==========================================================
        VkPresentInfoKHR presentInfo{};
        presentInfo.sType = VK_STRUCTURE_TYPE_PRESENT_INFO_KHR;
        presentInfo.waitSemaphoreCount = 1;
        presentInfo.pWaitSemaphores = signalSemaphores;

        VkSwapchainKHR swapchains[] = { m_swapchain->GetNativeSwapchain() };
        presentInfo.swapchainCount = 1;
        presentInfo.pSwapchains = swapchains;
        presentInfo.pImageIndices = &m_imageIndex;

        VkResult result = vkQueuePresentKHR(m_device->GetPresentQueue(), &presentInfo);

        // 如果呈现时发现窗口大小变了，重建 Swapchain
        if (result == VK_ERROR_OUT_OF_DATE_KHR || result == VK_SUBOPTIMAL_KHR) {
            RecreateSwapchain();
        } else if (result != VK_SUCCESS) {
            throw std::runtime_error("Failed to present swapchain image!");
        }

        m_isFrameStarted = false;
        m_currentFrame = (m_currentFrame + 1) % MAX_FRAMES_IN_FLIGHT;
    }

    void VulkanRenderer::BeginSwapChainRenderPass(VkCommandBuffer commandBuffer) {
        if (!m_isFrameStarted || commandBuffer != GetCurrentCommandBuffer()) {
            throw std::runtime_error("Can't begin render pass on an invalid command buffer!");
        }

        VkRenderPassBeginInfo renderPassInfo{};
        renderPassInfo.sType = VK_STRUCTURE_TYPE_RENDER_PASS_BEGIN_INFO;
        renderPassInfo.renderPass = m_defaultRenderPass->GetNativeRenderPass();
        renderPassInfo.framebuffer = m_defaultFramebuffers->GetNativeFramebuffers()[m_imageIndex];

        renderPassInfo.renderArea.offset = {0, 0};
        renderPassInfo.renderArea.extent = m_swapchain->GetExtent();

        // 黑灰色背景色
        std::array<VkClearValue, 2> clearValues{};
        clearValues[0].color = {{0.05f, 0.05f, 0.05f, 1.0f}};
        clearValues[1].depthStencil = {1.0f, 0}; // 如果以后加上深度测试，这里有用

        renderPassInfo.clearValueCount = static_cast<uint32_t>(clearValues.size());
        renderPassInfo.pClearValues = clearValues.data();

        vkCmdBeginRenderPass(commandBuffer, &renderPassInfo, VK_SUBPASS_CONTENTS_INLINE);

        // ★ 极其贴心的设计：底层自动帮你把 Viewport 和 Scissor 设好！业务层不用再操心！
        VkViewport viewport{};
        viewport.x = 0.0f;
        viewport.y = 0.0f;
        viewport.width = static_cast<float>(m_swapchain->GetExtent().width);
        viewport.height = static_cast<float>(m_swapchain->GetExtent().height);
        viewport.minDepth = 0.0f;
        viewport.maxDepth = 1.0f;
        VkRect2D scissor{{0, 0}, m_swapchain->GetExtent()};
        vkCmdSetViewport(commandBuffer, 0, 1, &viewport);
        vkCmdSetScissor(commandBuffer, 0, 1, &scissor);
    }

    void VulkanRenderer::EndSwapChainRenderPass(VkCommandBuffer commandBuffer) {
        if (!m_isFrameStarted || commandBuffer != GetCurrentCommandBuffer()) {
            throw std::runtime_error("Can't end render pass on an invalid command buffer!");
        }
        vkCmdEndRenderPass(commandBuffer);
    }

    VkCommandBuffer VulkanRenderer::GetCurrentCommandBuffer() const {
        if (!m_isFrameStarted) {
            throw std::runtime_error("Cannot get command buffer when frame not in progress!");
        }
        return m_commandBuffers[m_currentFrame]->GetNativeBuffer();
    }

} // namespace Lizeral
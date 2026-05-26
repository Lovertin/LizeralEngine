#pragma once

#include <string>
#include <vector>
#include <vulkan/vulkan.h>

namespace Lizeral {

    std::vector<char> ReadShaderFile(const std::string& filename);
    VkShaderModule CreateShaderModule(VkDevice device, const std::vector<char>& code);

} // namespace Lizeral

#version 460

layout(location = 0) in vec2 inUV;
layout(location = 0) out vec4 outColor; // 输出到当前的 History Buffer

// 我们需要 3 张输入贴图
layout(binding = 0) uniform sampler2D samplerSceneColor; // 当前帧 Lighting Pass 画出来的纯净画面
layout(binding = 1) uniform sampler2D samplerVelocity;   // 刚刚搞定的速度缓冲
layout(binding = 2) uniform sampler2D samplerHistory;    // 上一帧的 TAA 混合结果

void main() {
    // 1. 获取当前帧的颜色和速度
    vec4 currentColor = texture(samplerSceneColor, inUV);
    vec2 velocity = texture(samplerVelocity, inUV).xy;

    // 2. 重投影 (Reprojection)：计算这个像素在上一帧的 UV 坐标
    vec2 prevUV = inUV - velocity;

    // 如果越界（说明这个像素是这一帧刚从屏幕外进来的），直接抛弃历史，使用当前颜色
    if (prevUV.x < 0.0 || prevUV.x > 1.0 || prevUV.y < 0.0 || prevUV.y > 1.0) {
        outColor = currentColor;
        return;
    }

    // 3. 历史采样：去上一帧的画面里把颜色捞回来
    vec4 historyColor = texture(samplerHistory, prevUV);

    // ====================================================================
    // ★ 4. 邻域裁剪 (Neighborhood Clamping) - 绝对核心！没有它画面全是残影！
    // ====================================================================
    // 寻找当前像素周围 3x3 范围内的最亮和最暗颜色
    vec2 texelSize = 1.0 / textureSize(samplerSceneColor, 0);
    vec4 colorMin = currentColor;
    vec4 colorMax = currentColor;

    for (int x = -1; x <= 1; x++) {
        for (int y = -1; y <= 1; y++) {
            vec4 neighbor = texture(samplerSceneColor, inUV + vec2(x, y) * texelSize);
            colorMin = min(colorMin, neighbor);
            colorMax = max(colorMax, neighbor);
        }
    }

    // 强制把上一帧的颜色拉回到当前帧 3x3 邻域的色彩范围内
    // 如果不这么做，跑车开走后，留在原地的红色车身历史像素就会变成红色的“幽灵拖影”
    historyColor = clamp(historyColor, colorMin, colorMax);

    // ====================================================================
    // ★ 5. 时空混合 (Temporal Blend)
    // ====================================================================
    // 标准 TAA 混合比例：10% 相信当前帧，90% 相信历史帧（白嫖算力！）
    float blendWeight = 0.1; 
    
    // 当物体快速移动时（速度很大），我们可以增加对当前帧的信任，减少历史依赖，进一步防止拖影
    float velocityLength = length(velocity);
    if (velocityLength > 0.001) {
        blendWeight = mix(0.1, 0.4, clamp(velocityLength * 100.0, 0.0, 1.0));
    }

    outColor = mix(historyColor, currentColor, blendWeight);
}
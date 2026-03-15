#version 460
layout(location = 0) in vec2 inUV;

// 输出两张图！
layout(location = 0) out vec4 outGIHistory;    // RGB存GI颜色，Alpha存历史累加帧数 (History Length)
layout(location = 1) out vec2 outMoments;      // R存亮度一阶矩(均值)，G存二阶矩(平方均值)

layout(binding = 0) uniform sampler2D samplerNoisyGI;
layout(binding = 1) uniform sampler2D samplerNormalRoughness;
layout(binding = 2) uniform sampler2D samplerDepth;
layout(binding = 3) uniform sampler2D samplerVelocity;
layout(binding = 4) uniform sampler2D samplerPrevGIHistory;
layout(binding = 5) uniform sampler2D samplerPrevMoments;

// 提取亮度的辅助函数
float GetLuminance(vec3 color) {
    return dot(color, vec3(0.2126, 0.7152, 0.0722));
}

void main() {
    float depth = texture(samplerDepth, inUV).r;
    if (depth <= 0.0) {
        // 背景天空，直接输出当前值，清除历史
        outGIHistory = vec4(texture(samplerNoisyGI, inUV).rgb, 0.0);
        outMoments = vec2(0.0);
        return;
    }

    vec3 currentGI = texture(samplerNoisyGI, inUV).rgb;
    float currentLuma = GetLuminance(currentGI);
    
    vec2 velocity = texture(samplerVelocity, inUV).xy;
    vec2 prevUV = inUV - velocity;

    bool isHistoryValid = (prevUV.x >= 0.0 && prevUV.x <= 1.0 && prevUV.y >= 0.0 && prevUV.y <= 1.0);
    // TODO: 未来可以在这里加入上一帧的 Depth 和 Normal 校验，防止被遮挡物体的鬼影。目前先用 UV 校验。

    float historyLength = 0.0;
    vec3 prevGI = vec3(0.0);
    vec2 prevMoments = vec2(0.0);

    if (isHistoryValid) {
        vec4 prevGIHistory = texture(samplerPrevGIHistory, prevUV);
        prevGI = prevGIHistory.rgb;
        historyLength = prevGIHistory.a;
        prevMoments = texture(samplerPrevMoments, prevUV).rg;
    }
    
    vec3 giMin = currentGI;
    vec3 giMax = currentGI;
    vec2 texelSize = 1.0 / textureSize(samplerNoisyGI, 0);

    for (int y = -1; y <= 1; ++y) {
        for (int x = -1; x <= 1; ++x) {
            if (x == 0 && y == 0) continue;
            vec3 neighborGI = texture(samplerNoisyGI, inUV + vec2(x, y) * texelSize).rgb;
            giMin = min(giMin, neighborGI);
            giMax = max(giMax, neighborGI);
        }
    }

    // 将历史颜色强行限制在当前帧邻域的合法范围内！
    // 这样就算发生了视线遮挡，历史颜色也会瞬间被修正为当前帧的近似颜色。
    prevGI = clamp(prevGI, giMin, giMax);

    // 历史长度加1，最多累加 32 帧
    historyLength = min(32.0, isHistoryValid ? historyLength + 1.0 : 1.0);

    // 计算混合权重 (alpha)
    // 刚出现的物体 alpha=1(只相信当前帧)，稳定后的物体 alpha=0.05(更相信历史，极其平滑)
    float alpha = (historyLength < 4.0) ? 1.0 : max(0.05, 1.0 / historyLength);

    // 1. 时域混合 GI 颜色
    vec3 resolvedGI = mix(prevGI, currentGI, alpha);
    
    // 2. 时域混合 方差(Moments)
    vec2 currentMoments = vec2(currentLuma, currentLuma * currentLuma);
    vec2 resolvedMoments = mix(prevMoments, currentMoments, alpha);

    outGIHistory = vec4(resolvedGI, historyLength);
    outMoments = resolvedMoments;
}
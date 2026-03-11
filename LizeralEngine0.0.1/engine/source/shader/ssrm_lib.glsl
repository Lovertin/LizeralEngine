// ====================================================================
// Lizeral Engine - Screen Space Ray Marching 核心算法库
// ====================================================================

// 在这里，我们需要用到主 Shader 里绑定的 G-Buffer 采样器和 pc 矩阵
// 所以它们会由主 Shader 提供，我们直接使用

vec3 ScreenSpaceRayMarch(vec3 rayOrigin, vec3 rayDir, out float hitMask) {
    // 步进参数 (未来可以把这些暴露给 C++ 调节)
    const int   MAX_STEPS = 60;        // 最大步数
    const float STEP_SIZE = 0.2;       // 每次前进的 3D 距离
    const float THICKNESS = 0.5;       // 深度碰撞容差 (防止光线穿透物体打到背面的底面)

    vec3 currentPos = rayOrigin;
    hitMask = 0.0; // 0 表示没撞到，1 表示撞到了屏幕内的物体

    for(int i = 0; i < MAX_STEPS; i++) {
        // 1. 在 3D 世界中沿着射线往前走一步
        currentPos += rayDir * STEP_SIZE;

        // 2. 将当前的 3D 坐标“拍”回屏幕 2D 空间 (NDC)
        vec4 clipSpace = pc.viewProj * vec4(currentPos, 1.0);
        vec3 ndc = clipSpace.xyz / clipSpace.w;
        
        // NDC [-1, 1] 转 UV [0, 1]
        vec2 screenUV = ndc.xy * 0.5 + 0.5;

        // 3. 检查光线是否已经飞出了屏幕边界？飞出了就直接结束
        if(screenUV.x < 0.0 || screenUV.x > 1.0 || screenUV.y < 0.0 || screenUV.y > 1.0) {
            break;
        }

        // 4. 采样这一个像素在 G-Buffer 里的真实深度
        float sampledDepth = texture(samplerDepth, screenUV).r;

        // 5. 碰撞检测！
        // 在 Vulkan 中，NDC.z 也是 [0, 1]。
        // 如果我们当前光线的 Z (ndc.z) > 屏幕上记录的深度 (sampledDepth)，说明我们钻到物体“内部”或“背后”了！
        float depthDiff = ndc.z - sampledDepth;

        // 我们只接受深度差在一定厚度 (THICKNESS) 内的碰撞，避免误判打到远处的背景
        if(depthDiff > 0.0 && depthDiff < THICKNESS) {
            hitMask = 1.0; // 确认碰撞！
            
            // 撞到了！读取那个像素的 Albedo 颜色返回，这就是反射出来的颜色！
            return texture(samplerAlbedoMetallic, screenUV).rgb; 
        }
    }

    return vec3(0.0); // 什么都没撞到，返回纯黑
}
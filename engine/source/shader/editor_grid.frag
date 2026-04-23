#version 450

layout(push_constant) uniform PushConstants {
    mat4 invViewProj;
    vec4 cameraPosAndPlaneHeight;
    vec4 viewportSizeAndSpacing;
    vec4 fadeAndOpacity;
} pc;

layout(location = 0) out vec4 outColor;

vec3 homogeneousDivide(vec4 clipToWorld) {
    float w = clipToWorld.w;
    if (abs(w) < 1e-5) {
        w = (w < 0.0) ? -1e-5 : 1e-5;
    }
    return clipToWorld.xyz / w;
}

float gridLineFactor(vec2 worldPlane, float spacing) {
    vec2 scaled = worldPlane / max(spacing, 1e-5);
    vec2 derivatives = max(fwidth(scaled), vec2(1e-4));
    vec2 cell = abs(fract(scaled - 0.5) - 0.5) / derivatives;
    return 1.0 - min(min(cell.x, cell.y), 1.0);
}

float axisFactor(float coordinate) {
    float derivative = max(fwidth(coordinate), 1e-4);
    return 1.0 - min(abs(coordinate) / derivative, 1.0);
}

void main() {
    vec2 viewportSize = max(pc.viewportSizeAndSpacing.xy, vec2(1.0));
    vec2 uv = gl_FragCoord.xy / viewportSize;
    vec2 ndc = uv * 2.0 - 1.0;

    vec4 nearClip = vec4(ndc, 0.0, 1.0);
    vec4 farClip = vec4(ndc, 1.0, 1.0);

    vec4 nearWorld4 = pc.invViewProj * nearClip;
    vec4 farWorld4 = pc.invViewProj * farClip;
    vec3 nearWorld = homogeneousDivide(nearWorld4);
    vec3 farWorld = homogeneousDivide(farWorld4);

    vec3 rayDir = normalize(farWorld - nearWorld);
    float planeHeight = pc.cameraPosAndPlaneHeight.w;
    float denom = rayDir.y;

    if (abs(denom) < 1e-5) {
        discard;
    }

    float hitT = (planeHeight - nearWorld.y) / denom;
    if (hitT <= 0.0) {
        discard;
    }

    vec3 worldPos = nearWorld + rayDir * hitT;
    vec2 planePos = worldPos.xz;
    vec2 cameraPlanePos = pc.cameraPosAndPlaneHeight.xz;

    float distanceFade = 1.0 - clamp(length(planePos - cameraPlanePos) / pc.fadeAndOpacity.x, 0.0, 1.0);
    if (distanceFade <= 0.0) {
        discard;
    }

    float minorLine = gridLineFactor(planePos, pc.viewportSizeAndSpacing.z) * pc.fadeAndOpacity.z;
    float majorLine = gridLineFactor(planePos, pc.viewportSizeAndSpacing.w) * pc.fadeAndOpacity.y;
    float xAxisLine = axisFactor(worldPos.z) * pc.fadeAndOpacity.w;
    float zAxisLine = axisFactor(worldPos.x) * pc.fadeAndOpacity.w;

    vec3 color = vec3(0.0);
    color += vec3(0.35, 0.36, 0.39) * minorLine;
    color += vec3(0.62, 0.64, 0.68) * majorLine;
    color = mix(color, vec3(0.95, 0.22, 0.28), xAxisLine);
    color = mix(color, vec3(0.18, 0.56, 0.98), zAxisLine);

    float alpha = max(max(minorLine, majorLine), max(xAxisLine, zAxisLine)) * distanceFade;
    if (alpha <= 0.0) {
        discard;
    }

    outColor = vec4(color, clamp(alpha, 0.0, 1.0));
}

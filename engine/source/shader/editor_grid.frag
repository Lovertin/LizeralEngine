#version 450

layout(set = 0, binding = 0) uniform sampler2D uSceneDepth;

layout(push_constant) uniform PushConstants {
    mat4 viewProj;
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
    vec2 derivatives = max(fwidth(scaled), vec2(1e-3));
    vec2 cell = abs(fract(scaled - 0.5) - 0.5) / derivatives;
    float line = 1.0 - min(min(cell.x, cell.y), 1.0);
    return smoothstep(0.0, 1.0, line);
}

float axisFactor(float coordinate, float minorSpacing) {
    float derivative = max(fwidth(coordinate), minorSpacing * 0.02);
    return 1.0 - smoothstep(0.0, derivative, abs(coordinate));
}

float projectDepth(vec3 worldPos) {
    vec4 clip = pc.viewProj * vec4(worldPos, 1.0);
    float w = clip.w;
    if (abs(w) < 1e-5) {
        w = (w < 0.0) ? -1e-5 : 1e-5;
    }
    return clip.z / w;
}

void main() {
    vec2 viewportSize = max(pc.viewportSizeAndSpacing.xy, vec2(1.0));
    vec2 uv = gl_FragCoord.xy / viewportSize;
    vec2 ndc = uv * 2.0 - 1.0;

    vec4 farClip = vec4(ndc, 0.0, 1.0);

    vec4 farWorld4 = pc.invViewProj * farClip;
    vec3 cameraPos = pc.cameraPosAndPlaneHeight.xyz;
    vec3 farWorld = homogeneousDivide(farWorld4);

    vec3 rayDir = normalize(farWorld - cameraPos);
    float planeHeight = pc.cameraPosAndPlaneHeight.w;
    float denom = rayDir.y;
    bool hasGridHit = abs(denom) >= 1e-5;
    float hitT = hasGridHit ? (planeHeight - cameraPos.y) / denom : -1.0;
    hasGridHit = hasGridHit && hitT > 0.0;

    vec3 worldPos = vec3(0.0);
    vec2 planePos = vec2(0.0);
    vec2 cameraPlanePos = cameraPos.xz;
    float distanceFade = 0.0;

    if (hasGridHit) {
        worldPos = cameraPos + rayDir * hitT;
        planePos = worldPos.xz;
        distanceFade = 1.0 - clamp(length(planePos - cameraPlanePos) / pc.fadeAndOpacity.x, 0.0, 1.0);
    }

    float minorLine = 0.0;
    float majorLine = 0.0;
    if (hasGridHit && distanceFade > 0.0) {
        minorLine = gridLineFactor(planePos, pc.viewportSizeAndSpacing.z) * pc.fadeAndOpacity.z;
        majorLine = gridLineFactor(planePos, pc.viewportSizeAndSpacing.w) * pc.fadeAndOpacity.y;
    }

    float xAxisLine = 0.0;
    float zAxisLine = 0.0;
    if (hasGridHit) {
        xAxisLine = axisFactor(worldPos.z, pc.viewportSizeAndSpacing.z) * pc.fadeAndOpacity.w;
        zAxisLine = axisFactor(worldPos.x, pc.viewportSizeAndSpacing.z) * pc.fadeAndOpacity.w;
    }

    vec3 color = vec3(0.0);
    color += vec3(0.35, 0.36, 0.39) * minorLine;
    color += vec3(0.62, 0.64, 0.68) * majorLine;
    color = mix(color, vec3(0.95, 0.22, 0.28), xAxisLine);
    color = mix(color, vec3(0.18, 0.56, 0.98), zAxisLine);

    float alpha = max(max(minorLine, majorLine), max(xAxisLine, zAxisLine)) * max(distanceFade, 0.0);
    if (alpha <= 0.0) {
        discard;
    }

    float sceneDepth = texture(uSceneDepth, uv).r;
    float gridDepth = projectDepth(worldPos);
    if (sceneDepth > gridDepth + 1e-4) {
        discard;
    }

    outColor = vec4(color, clamp(alpha, 0.0, 1.0));
}

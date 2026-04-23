#version 450

layout(set = 0, binding = 0) uniform sampler2D uSceneDepth;

layout(push_constant) uniform PushConstants {
    mat4 viewProj;
    vec4 viewportSizeAndThickness;
    vec4 axisStart;
    vec4 axisEnd;
    vec4 axisColor;
} pc;

layout(location = 0) out vec4 outColor;
layout(depth_greater) out float gl_FragDepth;

const float kMinClipW = 1e-4;

bool clipAgainstPositivePlane(float distanceA, float distanceB, inout vec4 pointA, inout vec4 pointB) {
    bool aInside = distanceA >= 0.0;
    bool bInside = distanceB >= 0.0;
    if (!aInside && !bInside) {
        return false;
    }
    if (aInside && bInside) {
        return true;
    }

    float denom = distanceA - distanceB;
    float t = (abs(denom) > 1e-6) ? clamp(distanceA / denom, 0.0, 1.0) : 0.0;
    vec4 intersection = mix(pointA, pointB, t);
    if (!aInside) {
        pointA = intersection;
    } else {
        pointB = intersection;
    }
    return true;
}

bool clipAxisSegment(inout vec4 clipA, inout vec4 clipB) {
    if (!clipAgainstPositivePlane(clipA.w - kMinClipW, clipB.w - kMinClipW, clipA, clipB)) {
        return false;
    }
    if (!clipAgainstPositivePlane(clipA.z, clipB.z, clipA, clipB)) {
        return false;
    }
    if (!clipAgainstPositivePlane(clipA.w - clipA.z, clipB.w - clipB.z, clipA, clipB)) {
        return false;
    }
    return true;
}

vec2 projectScreen(vec4 clipPos, out float depth) {
    float safeW = abs(clipPos.w) < kMinClipW ? ((clipPos.w < 0.0) ? -kMinClipW : kMinClipW) : clipPos.w;
    vec2 ndc = clipPos.xy / safeW;
    depth = clipPos.z / safeW;
    return (ndc * 0.5 + 0.5) * pc.viewportSizeAndThickness.xy;
}

void main() {
    vec2 viewportSize = max(pc.viewportSizeAndThickness.xy, vec2(1.0));
    vec2 uv = gl_FragCoord.xy / viewportSize;

    vec4 clipA = pc.axisStart;
    vec4 clipB = pc.axisEnd;
    if (!clipAxisSegment(clipA, clipB)) {
        discard;
    }

    float depthA = 0.0;
    float depthB = 0.0;
    vec2 a = projectScreen(clipA, depthA);
    vec2 b = projectScreen(clipB, depthB);

    vec2 ab = b - a;
    float abLenSq = max(dot(ab, ab), 1e-5);
    float t = clamp(dot(gl_FragCoord.xy - a, ab) / abLenSq, 0.0, 1.0);
    vec2 closest = mix(a, b, t);
    float distancePx = length(gl_FragCoord.xy - closest);

    float thicknessPx = max(pc.viewportSizeAndThickness.z, 1.0);
    float alpha = 1.0 - smoothstep(thicknessPx, thicknessPx + 1.25, distancePx);
    if (alpha <= 0.0) {
        discard;
    }

    float axisDepth = mix(depthA, depthB, t);
    float sceneDepth = texture(uSceneDepth, uv).r;
    if (sceneDepth > axisDepth + 1e-4) {
        discard;
    }

    gl_FragDepth = clamp(axisDepth, 0.0, 1.0);
    outColor = vec4(pc.axisColor.rgb, alpha * pc.axisColor.a);
}

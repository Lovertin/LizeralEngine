#version 450

// DebugLineVertex position (location = 0)
layout(location = 0) in vec3 inPosition;
// DebugLineVertex color (location = 1)
layout(location = 1) in vec3 inColor;

// Push Constants

layout(push_constant) uniform PushConstants {
    mat4 viewProj;
} pc;

layout(location = 0) out vec3 fragColor;
layout(location = 1) out vec2 fragUv;
layout(location = 2) out float fragDepth;

void main() {
    gl_Position = pc.viewProj * vec4(inPosition, 1.0); //unjetter

    fragColor = inColor;

    float safeW = gl_Position.w;
    if (abs(safeW) < 1e-5) {
        safeW = (safeW < 0.0) ? -1e-5 : 1e-5;
    }
    vec2 ndc = gl_Position.xy / safeW;
    fragUv = ndc * 0.5 + 0.5;
    fragDepth = gl_Position.z / safeW;
}

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

void main() {
    gl_Position = pc.viewProj * vec4(inPosition, 1.0); //unjetter

    fragColor = inColor;
}
#version 450

layout(set = 0, binding = 0) uniform sampler2D uSceneDepth;

layout(location = 0) in vec3 fragColor;
layout(location = 1) in vec2 fragUv;
layout(location = 2) in float fragDepth;

layout(location = 0) out vec4 outColor;

void main() {
    float sceneDepth = texture(uSceneDepth, fragUv).r;
    if (sceneDepth > fragDepth + 1e-4) {
        discard;
    }

    outColor = vec4(fragColor, 1.0);
}

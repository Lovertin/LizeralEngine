#version 460
layout(location = 0) in vec2 inUV;
layout(location = 0) out vec4 outColor;

layout(binding = 0) uniform sampler2D samplerColor;

vec3 ACESFilm(vec3 x) {
    float a = 2.51f;
    float b = 0.03f;
    float c = 2.43f;
    float d = 0.59f;
    float e = 0.14f;
    return clamp((x*(a*x+b))/(x*(c*x+d)+e), 0.0, 1.0);
}

void main() {

    vec3 hdrColor = texture(samplerColor, inUV).rgb;
    
    float exposure = 1.0; 
    hdrColor *= exposure;

    vec3 mappedColor = ACESFilm(hdrColor);

    mappedColor = pow(mappedColor, vec3(1.0 / 2.2));

    outColor = vec4(mappedColor, 1.0);
}
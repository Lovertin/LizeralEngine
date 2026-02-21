#version 330 core
out vec4 FragColor;

in vec3 TexCoords;

uniform samplerCube skybox; // 输入 Cubemap

void main()
{    
    // 直接采样 Cubemap。因为是 HDR 图，可能需要简单的色调映射防止过曝
    vec3 envColor = texture(skybox, TexCoords).rgb;
    
    // 简单的 Reinhard 色调映射 + Gamma 校正 (让 HDR 显示到屏幕上更自然)
    // 如果你觉得天空太暗或太亮，可以调整这里的逻辑
    envColor = envColor / (envColor + vec3(1.0));
    envColor = pow(envColor, vec3(1.0/2.2)); 
  
    FragColor = vec4(envColor, 1.0);
}
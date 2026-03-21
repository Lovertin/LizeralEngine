
// Lizeral Engine - Screen Space Ray Marching


vec3 ScreenSpaceRayMarch(vec3 rayOrigin, vec3 rayDir, out float hitMask) {
    // TODO: expose to cpp
    const int   MAX_STEPS = 60;        // Max Step
    const float STEP_SIZE = 0.2;       
    const float THICKNESS = 0.5;       // bias of each Collision

    vec3 currentPos = rayOrigin;
    hitMask = 0.0; 

    for(int i = 0; i < MAX_STEPS; i++) {
        currentPos += rayDir * STEP_SIZE;

        vec4 clipSpace = pc.viewProj * vec4(currentPos, 1.0);
        vec3 ndc = clipSpace.xyz / clipSpace.w;

        vec2 screenUV = ndc.xy * 0.5 + 0.5;

        if(screenUV.x < 0.0 || screenUV.x > 1.0 || screenUV.y < 0.0 || screenUV.y > 1.0) {
            break;
        }

        float sampledDepth = texture(samplerDepth, screenUV).r;

        float depthDiff = ndc.z - sampledDepth;

        // retain the Diff less than THICKNESS
        if(depthDiff > 0.0 && depthDiff < THICKNESS) {
            hitMask = 1.0; // catch it !!
            return texture(samplerAlbedoMetallic, screenUV).rgb; 
        }
    }

    return vec3(0.0); // return null
}
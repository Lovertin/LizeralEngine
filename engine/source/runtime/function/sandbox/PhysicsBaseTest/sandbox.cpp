#include <cmath>
#include <exception>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#define GLFW_INCLUDE_VULKAN
#include <GLFW/glfw3.h>

#include "runtime/function/render/rhi/vulkan/VulkanContext.h"
#include "runtime/function/render/rhi/vulkan/VulkanDevice.h"
#include "runtime/function/render/VulkanRenderer/VulkanRenderer.h"
#include "runtime/function/render/VulkanRenderingSystem/VulkanRenderingSystem.h"

#include "runtime/function/ecs/registry.h"
#include "runtime/function/ecs/hybrid/hybrid_registry.h"
#include "runtime/function/ecs/hybrid/hybrid_default_traits.h"
#include "runtime/function/ecs/components/Transform/TransformComponent.h"
#include "runtime/function/ecs/components/Camera/CameraComponent.h"
#include "runtime/function/ecs/components/Camera/CameraControlComponent.h"
#include "runtime/function/ecs/components/Collider/ColliderComponent.h"
#include "runtime/function/ecs/components/RigidBody/RigidBodyComponent.h"
#include "runtime/function/input/input.h"
#include "runtime/function/render/CameraControlSystem/CameraControlSystem.h"
#include "runtime/function/render/CameraSystem/CameraSystem.h"
#include "runtime/function/physics/PhysicsScene.h"
#include "runtime/function/physics/PhysicsSystem.h"

using namespace Lizeral;

namespace {

    constexpr uint32_t kWindowWidth = 1280;
    constexpr uint32_t kWindowHeight = 720;

    constexpr int kGridX = 10;
    constexpr int kGridZ = 10;
    constexpr int kStackHeight = 6;
    constexpr int kMoverCount = 12;
    constexpr int kWarmupFrames = 120;
    constexpr int kSampleWindow = 120;

    struct MirrorEntitySet {
        Entity sparse = null_entity;
        Entity hybridView = null_entity;
        Entity hybridChunk = null_entity;
    };

    struct KinematicMoverPair {
        MirrorEntitySet entities;
        Vector3 basePosition;
        float amplitude = 0.0f;
        float speed = 0.0f;
        float phase = 0.0f;
    };

    struct BenchmarkWindow {
        int sampleCount = 0;
        PhysicsAccessProfile sparseAccess{};
        PhysicsAccessProfile hybridViewAccess{};
        PhysicsAccessProfile hybridChunkAccess{};
        PhysicsTickProfile sparseTick{};
        PhysicsTickProfile hybridViewTick{};
        PhysicsTickProfile hybridChunkTick{};

        void Add(
            const PhysicsAccessProfile& sparseAccessProfile,
            const PhysicsAccessProfile& hybridViewAccessProfile,
            const PhysicsAccessProfile& hybridChunkAccessProfile,
            const PhysicsTickProfile& sparse,
            const PhysicsTickProfile& hybridView,
            const PhysicsTickProfile& hybridChunk)
        {
            ++sampleCount;

            sparseAccess.entity_count = sparseAccessProfile.entity_count;
            sparseAccess.view_build_us += sparseAccessProfile.view_build_us;
            sparseAccess.iteration_us += sparseAccessProfile.iteration_us;
            sparseAccess.total_us += sparseAccessProfile.total_us;

            hybridViewAccess.entity_count = hybridViewAccessProfile.entity_count;
            hybridViewAccess.view_build_us += hybridViewAccessProfile.view_build_us;
            hybridViewAccess.iteration_us += hybridViewAccessProfile.iteration_us;
            hybridViewAccess.total_us += hybridViewAccessProfile.total_us;

            hybridChunkAccess.entity_count = hybridChunkAccessProfile.entity_count;
            hybridChunkAccess.view_build_us += hybridChunkAccessProfile.view_build_us;
            hybridChunkAccess.iteration_us += hybridChunkAccessProfile.iteration_us;
            hybridChunkAccess.total_us += hybridChunkAccessProfile.total_us;

            sparseTick.entity_count = sparse.entity_count;
            sparseTick.created_body_count += sparse.created_body_count;
            sparseTick.dirty_sync_count += sparse.dirty_sync_count;
            sparseTick.transform_write_back_count += sparse.transform_write_back_count;
            sparseTick.pre_step_us += sparse.pre_step_us;
            sparseTick.simulation_us += sparse.simulation_us;
            sparseTick.write_back_us += sparse.write_back_us;
            sparseTick.total_us += sparse.total_us;

            hybridViewTick.entity_count = hybridView.entity_count;
            hybridViewTick.created_body_count += hybridView.created_body_count;
            hybridViewTick.dirty_sync_count += hybridView.dirty_sync_count;
            hybridViewTick.transform_write_back_count += hybridView.transform_write_back_count;
            hybridViewTick.pre_step_us += hybridView.pre_step_us;
            hybridViewTick.simulation_us += hybridView.simulation_us;
            hybridViewTick.write_back_us += hybridView.write_back_us;
            hybridViewTick.total_us += hybridView.total_us;

            hybridChunkTick.entity_count = hybridChunk.entity_count;
            hybridChunkTick.created_body_count += hybridChunk.created_body_count;
            hybridChunkTick.dirty_sync_count += hybridChunk.dirty_sync_count;
            hybridChunkTick.transform_write_back_count += hybridChunk.transform_write_back_count;
            hybridChunkTick.pre_step_us += hybridChunk.pre_step_us;
            hybridChunkTick.simulation_us += hybridChunk.simulation_us;
            hybridChunkTick.write_back_us += hybridChunk.write_back_us;
            hybridChunkTick.total_us += hybridChunk.total_us;
        }

        void Reset() {
            sampleCount = 0;
            sparseAccess = {};
            hybridViewAccess = {};
            hybridChunkAccess = {};
            sparseTick = {};
            hybridViewTick = {};
            hybridChunkTick = {};
        }
    };

    double Average(double total, int count) {
        return count > 0 ? total / static_cast<double>(count) : 0.0;
    }

    double Speedup(double baseline, double candidate) {
        return candidate > 0.0 ? baseline / candidate : 0.0;
    }

    template<typename RegistryT>
    Entity SpawnPhysicsBox(
        RegistryT& registry,
        const Vector3& position,
        const Vector3& colliderSize,
        float mass,
        bool kinematic,
        bool showDebug,
        const Vector3& scale = Vector3::UNIT_SCALE)
    {
        const Entity entity = registry.create();

        auto& transform = registry.template emplace<TransformComponent>(entity);
        transform.setPosition(position);
        transform.setScale(scale);

        auto& rigidBody = registry.template emplace<RigidBodyComponent>(entity);
        rigidBody.setMass(mass);
        rigidBody.setKinematic(kinematic);
        rigidBody.setFriction(0.8f);
        rigidBody.setRestitution(0.05f);

        auto& collider = registry.template emplace<ColliderComponent>(entity);
        collider.setType(ColliderType::Box);
        collider.setSize(colliderSize);
        collider.setDebug(showDebug);

        return entity;
    }

    void SpawnMirrorPhysicsBox(
        Registry& sparseRegistry,
        HybridRegistry& hybridViewRegistry,
        HybridRegistry& hybridChunkRegistry,
        const Vector3& position,
        const Vector3& colliderSize,
        float mass,
        bool kinematic,
        bool showDebug)
    {
        SpawnPhysicsBox(sparseRegistry, position, colliderSize, mass, kinematic, showDebug);
        SpawnPhysicsBox(hybridViewRegistry, position, colliderSize, mass, kinematic, showDebug);
        SpawnPhysicsBox(hybridChunkRegistry, position, colliderSize, mass, kinematic, showDebug);
    }

    void BuildBenchmarkScene(
        Registry& sparseRegistry,
        HybridRegistry& hybridViewRegistry,
        HybridRegistry& hybridChunkRegistry,
        std::vector<KinematicMoverPair>& movers)
    {
        SpawnMirrorPhysicsBox(
            sparseRegistry,
            hybridViewRegistry,
            hybridChunkRegistry,
            Vector3(0.0f, -0.75f, 0.0f),
            Vector3(80.0f, 1.5f, 80.0f),
            0.0f,
            false,
            true);

        const Vector3 boxSize(0.9f, 0.9f, 0.9f);
        const float spacing = 1.1f;
        const float startX = -0.5f * static_cast<float>(kGridX - 1) * spacing;
        const float startZ = -0.5f * static_cast<float>(kGridZ - 1) * spacing;

        for (int x = 0; x < kGridX; ++x) {
            for (int z = 0; z < kGridZ; ++z) {
                for (int y = 0; y < kStackHeight; ++y) {
                    const Vector3 position(
                        startX + static_cast<float>(x) * spacing,
                        0.6f + static_cast<float>(y) * spacing,
                        startZ + static_cast<float>(z) * spacing);

                    SpawnMirrorPhysicsBox(
                        sparseRegistry,
                        hybridViewRegistry,
                        hybridChunkRegistry,
                        position,
                        boxSize,
                        1.0f,
                        false,
                        true);
                }
            }
        }

        movers.reserve(kMoverCount);
        for (int i = 0; i < kMoverCount; ++i) {
            const float normalized = static_cast<float>(i) / static_cast<float>(kMoverCount);
            const Vector3 basePosition(-18.0f + normalized * 36.0f, 1.5f, -14.0f);

            KinematicMoverPair mover{};
            mover.basePosition = basePosition;
            mover.amplitude = 2.5f + 0.15f * static_cast<float>(i);
            mover.speed = 1.0f + 0.08f * static_cast<float>(i);
            mover.phase = normalized * 3.1415926f;

            mover.entities.sparse = SpawnPhysicsBox(
                sparseRegistry,
                basePosition,
                Vector3(2.2f, 0.5f, 2.2f),
                0.0f,
                true,
                true);

            mover.entities.hybridView = SpawnPhysicsBox(
                hybridViewRegistry,
                basePosition,
                Vector3(2.2f, 0.5f, 2.2f),
                0.0f,
                true,
                true);

            mover.entities.hybridChunk = SpawnPhysicsBox(
                hybridChunkRegistry,
                basePosition,
                Vector3(2.2f, 0.5f, 2.2f),
                0.0f,
                true,
                true);

            movers.push_back(mover);
        }
    }

    void CreateSparseCamera(Registry& registry) {
        const Entity cameraEntity = registry.create();
        auto& transform = registry.emplace<TransformComponent>(cameraEntity);
        transform.setPosition(Vector3(0.0f, 14.0f, 28.0f));
        transform.setForward(Vector3(0.0f, -0.3f, -1.0f));

        auto& camera = registry.emplace<CameraComponent>(cameraEntity);
        camera.setFov(45.0f);
        camera.setAspect(static_cast<float>(kWindowWidth) / static_cast<float>(kWindowHeight));
        camera.setzNear(0.001f);
        camera.setzFar(1000.0f);
        camera.setMain(true);

        auto& control = registry.emplace<CameraControlComponent>(cameraEntity);
        control.setMoveSpeed(12.0f);
        control.setSensitivityX(0.1f);
        control.setSensitivityY(0.1f);
        control.setSpeedMutiplier(5.0f);
    }

    void UpdateKinematicMovers(
        float timeSeconds,
        Registry& sparseRegistry,
        HybridRegistry& hybridViewRegistry,
        HybridRegistry& hybridChunkRegistry,
        const std::vector<KinematicMoverPair>& movers)
    {
        for (const KinematicMoverPair& mover : movers) {
            const float offsetX = std::sin(timeSeconds * mover.speed + mover.phase) * mover.amplitude;
            const float offsetY = std::cos(timeSeconds * (mover.speed * 0.5f) + mover.phase) * 0.35f;
            const Vector3 newPosition = mover.basePosition + Vector3(offsetX, offsetY, 0.0f);

            sparseRegistry.get<TransformComponent>(mover.entities.sparse).setPosition(newPosition);
            hybridViewRegistry.get<TransformComponent>(mover.entities.hybridView).setPosition(newPosition);
            hybridChunkRegistry.get<TransformComponent>(mover.entities.hybridChunk).setPosition(newPosition);
        }
    }

    void PrintBenchmarkWindow(const BenchmarkWindow& window) {
        if (window.sampleCount <= 0) {
            return;
        }

        const double sparseAccessTotal = Average(window.sparseAccess.total_us, window.sampleCount);
        const double hybridViewAccessTotal = Average(window.hybridViewAccess.total_us, window.sampleCount);
        const double hybridChunkAccessTotal = Average(window.hybridChunkAccess.total_us, window.sampleCount);
        const double sparseAccessView = Average(window.sparseAccess.view_build_us, window.sampleCount);
        const double hybridViewAccessView = Average(window.hybridViewAccess.view_build_us, window.sampleCount);
        const double hybridChunkAccessView = Average(window.hybridChunkAccess.view_build_us, window.sampleCount);
        const double sparseAccessIter = Average(window.sparseAccess.iteration_us, window.sampleCount);
        const double hybridViewAccessIter = Average(window.hybridViewAccess.iteration_us, window.sampleCount);
        const double hybridChunkAccessIter = Average(window.hybridChunkAccess.iteration_us, window.sampleCount);

        const double sparseTickPre = Average(window.sparseTick.pre_step_us, window.sampleCount);
        const double sparseTickSim = Average(window.sparseTick.simulation_us, window.sampleCount);
        const double sparseTickWrite = Average(window.sparseTick.write_back_us, window.sampleCount);
        const double sparseTickTotal = Average(window.sparseTick.total_us, window.sampleCount);

        const double hybridViewTickPre = Average(window.hybridViewTick.pre_step_us, window.sampleCount);
        const double hybridViewTickSim = Average(window.hybridViewTick.simulation_us, window.sampleCount);
        const double hybridViewTickWrite = Average(window.hybridViewTick.write_back_us, window.sampleCount);
        const double hybridViewTickTotal = Average(window.hybridViewTick.total_us, window.sampleCount);

        const double hybridChunkTickPre = Average(window.hybridChunkTick.pre_step_us, window.sampleCount);
        const double hybridChunkTickSim = Average(window.hybridChunkTick.simulation_us, window.sampleCount);
        const double hybridChunkTickWrite = Average(window.hybridChunkTick.write_back_us, window.sampleCount);
        const double hybridChunkTickTotal = Average(window.hybridChunkTick.total_us, window.sampleCount);

        std::cout << "\n[PhysicsBenchmark] Samples=" << window.sampleCount
                  << " Entities=" << window.sparseTick.entity_count << '\n';

        std::cout << std::fixed << std::setprecision(2);
        std::cout << "  Access Sparse  : total=" << sparseAccessTotal
                  << "us view=" << sparseAccessView
                  << "us iter=" << sparseAccessIter << "us\n";
        std::cout << "  Access H-View  : total=" << hybridViewAccessTotal
                  << "us view=" << hybridViewAccessView
                  << "us iter=" << hybridViewAccessIter
                  << "us speedup=" << Speedup(sparseAccessTotal, hybridViewAccessTotal) << "x\n";
        std::cout << "  Access H-Chunk : total=" << hybridChunkAccessTotal
                  << "us view=" << hybridChunkAccessView
                  << "us iter=" << hybridChunkAccessIter
                  << "us sparse-speedup=" << Speedup(sparseAccessTotal, hybridChunkAccessTotal)
                  << "x view-speedup=" << Speedup(hybridViewAccessTotal, hybridChunkAccessTotal) << "x\n";

        std::cout << "  Tick Sparse    : pre=" << sparseTickPre
                  << "us sim=" << sparseTickSim
                  << "us write=" << sparseTickWrite
                  << "us total=" << sparseTickTotal << "us\n";
        std::cout << "  Tick H-View    : pre=" << hybridViewTickPre
                  << "us sim=" << hybridViewTickSim
                  << "us write=" << hybridViewTickWrite
                  << "us total=" << hybridViewTickTotal
                  << "us total-speedup=" << Speedup(sparseTickTotal, hybridViewTickTotal) << "x\n";
        std::cout << "  Tick H-Chunk   : pre=" << hybridChunkTickPre
                  << "us sim=" << hybridChunkTickSim
                  << "us write=" << hybridChunkTickWrite
                  << "us total=" << hybridChunkTickTotal
                  << "us sparse-speedup=" << Speedup(sparseTickTotal, hybridChunkTickTotal)
                  << "x view-speedup=" << Speedup(hybridViewTickTotal, hybridChunkTickTotal) << "x\n";

        std::cout << "  Sparse Counters: created=" << window.sparseTick.created_body_count
                  << " dirty=" << window.sparseTick.dirty_sync_count
                  << " writeback=" << window.sparseTick.transform_write_back_count << '\n';
        std::cout << "  H-View Counters: created=" << window.hybridViewTick.created_body_count
                  << " dirty=" << window.hybridViewTick.dirty_sync_count
                  << " writeback=" << window.hybridViewTick.transform_write_back_count << '\n';
        std::cout << "  H-Chunk Count. : created=" << window.hybridChunkTick.created_body_count
                  << " dirty=" << window.hybridChunkTick.dirty_sync_count
                  << " writeback=" << window.hybridChunkTick.transform_write_back_count << '\n';
    }

    void glfwKeyCallback(GLFWwindow*, int key, int, int action, int) {
        Key lizKey;
        bool found = true;
        switch (key) {
            case GLFW_KEY_W: lizKey = Key::W; break;
            case GLFW_KEY_A: lizKey = Key::A; break;
            case GLFW_KEY_S: lizKey = Key::S; break;
            case GLFW_KEY_D: lizKey = Key::D; break;
            case GLFW_KEY_Q: lizKey = Key::Q; break;
            case GLFW_KEY_E: lizKey = Key::E; break;
            case GLFW_KEY_LEFT_SHIFT: lizKey = Key::LEFT_SHIFT; break;
            case GLFW_KEY_ESCAPE: lizKey = Key::ESC; break;
            default: found = false; break;
        }

        if (found) {
            Input::GetInstance().SetKeyDown(lizKey, action != GLFW_RELEASE);
        }
    }

    void glfwCursorPosCallback(GLFWwindow*, double xpos, double ypos) {
        Input::GetInstance().SetMousePosition(static_cast<float>(xpos), static_cast<float>(ypos));
    }

    void glfwMouseButtonCallback(GLFWwindow*, int button, int action, int) {
        if (button == GLFW_MOUSE_BUTTON_RIGHT) {
            Input::GetInstance().SetMouseButtonDown(MouseButton::Right, action != GLFW_RELEASE);
        }
    }

} // namespace

int main() {
    GLFWwindow* window = nullptr;
    VkSurfaceKHR surface = VK_NULL_HANDLE;
    VkInstance instance = VK_NULL_HANDLE;

    try {
        glfwInit();
        glfwWindowHint(GLFW_CLIENT_API, GLFW_NO_API);
        glfwWindowHint(GLFW_RESIZABLE, GLFW_TRUE);
        window = glfwCreateWindow(
            kWindowWidth,
            kWindowHeight,
            "Lizeral Physics Benchmark - SparseSet vs Hybrid Archetype",
            nullptr,
            nullptr);

        glfwSetKeyCallback(window, glfwKeyCallback);
        glfwSetCursorPosCallback(window, glfwCursorPosCallback);
        glfwSetMouseButtonCallback(window, glfwMouseButtonCallback);

        VulkanContext vulkanContext;
        uint32_t glfwExtensionCount = 0;
        const char** glfwExtensions = glfwGetRequiredInstanceExtensions(&glfwExtensionCount);
        std::vector<const char*> requiredExtensions(glfwExtensions, glfwExtensions + glfwExtensionCount);
        vulkanContext.Initialize("Lizeral Physics Benchmark", requiredExtensions);
        instance = (VkInstance)vulkanContext.GetNativeInstance();

        if (glfwCreateWindowSurface(instance, window, nullptr, &surface) != VK_SUCCESS) {
            throw std::runtime_error("Failed to create window surface.");
        }

        {
            VulkanDevice vulkanDevice(&vulkanContext, surface);
            VulkanRenderer renderer(&vulkanContext, &vulkanDevice, kWindowWidth, kWindowHeight);
            VulkanRenderingSystem renderSystem;

            Registry sparseRegistry;
            HybridRegistry hybridViewRegistry;
            HybridRegistry hybridChunkRegistry;
            PhysicsScene sparseScene;
            PhysicsScene hybridViewScene;
            PhysicsScene hybridChunkScene;
            PhysicsSystem sparsePhysics;
            PhysicsSystem hybridViewPhysics;
            PhysicsSystem hybridChunkPhysics;
            CameraControlSystem cameraControlSystem;
            CameraSystem cameraSystem;
            std::vector<KinematicMoverPair> movers;
            BenchmarkWindow benchmarkWindow;

            sparseScene.Initialize();
            hybridViewScene.Initialize();
            hybridChunkScene.Initialize();
            sparsePhysics.Initialize(&sparseScene);
            hybridViewPhysics.Initialize(&hybridViewScene);
            hybridChunkPhysics.Initialize(&hybridChunkScene);
            hybridViewPhysics.SetHybridExecutionMode(PhysicsHybridExecutionMode::View);
            hybridChunkPhysics.SetHybridExecutionMode(PhysicsHybridExecutionMode::Chunked);

            CreateSparseCamera(sparseRegistry);
            BuildBenchmarkScene(sparseRegistry, hybridViewRegistry, hybridChunkRegistry, movers);

            VkExtent2D actualExtent = renderer.GetSwapchainExtent();
            renderSystem.Initialize(&vulkanContext, &vulkanDevice, &renderer, actualExtent.width, actualExtent.height);
            renderSystem.SetViewport(0, 0, actualExtent.width, actualExtent.height);

            std::cout << "\n[PhysicsBenchmark] Warmup frames: " << kWarmupFrames
                      << " | Sample window: " << kSampleWindow
                      << " | Dynamic bodies: " << (kGridX * kGridZ * kStackHeight)
                      << " | Movers: " << kMoverCount << std::endl;
            std::cout << "[PhysicsBenchmark] Hold RMB + WASD to inspect the sparse-set world." << std::endl;

            float lastTime = static_cast<float>(glfwGetTime());
            int frameIndex = 0;

            while (!glfwWindowShouldClose(window)) {
                const float currentTime = static_cast<float>(glfwGetTime());
                const float dt = currentTime - lastTime;
                lastTime = currentTime;

                Input::GetInstance().Tick();
                if (Input::GetInstance().GetKey(Key::ESC)) {
                    glfwSetWindowShouldClose(window, true);
                }
                glfwPollEvents();

                cameraControlSystem.Tick(dt, sparseRegistry);
                cameraSystem.Tick(sparseRegistry);

                UpdateKinematicMovers(currentTime, sparseRegistry, hybridViewRegistry, hybridChunkRegistry, movers);

                if (frameIndex < kWarmupFrames) {
                    sparsePhysics.Tick(dt, sparseRegistry);
                    hybridViewPhysics.Tick(dt, hybridViewRegistry);
                    hybridChunkPhysics.Tick(dt, hybridChunkRegistry);
                } else {
                    const PhysicsAccessProfile sparseAccess = sparsePhysics.ProfileComponentAccess(sparseRegistry);
                    const PhysicsAccessProfile hybridViewAccess = hybridViewPhysics.ProfileComponentAccess(hybridViewRegistry);
                    const PhysicsAccessProfile hybridChunkAccess = hybridChunkPhysics.ProfileChunkAccess(hybridChunkRegistry);
                    const PhysicsTickProfile sparseTick = sparsePhysics.TickProfiled(dt, sparseRegistry);
                    const PhysicsTickProfile hybridViewTick = hybridViewPhysics.TickProfiled(dt, hybridViewRegistry);
                    const PhysicsTickProfile hybridChunkTick = hybridChunkPhysics.TickProfiled(dt, hybridChunkRegistry);

                    benchmarkWindow.Add(
                        sparseAccess,
                        hybridViewAccess,
                        hybridChunkAccess,
                        sparseTick,
                        hybridViewTick,
                        hybridChunkTick);

                    if (benchmarkWindow.sampleCount >= kSampleWindow) {
                        PrintBenchmarkWindow(benchmarkWindow);
                        benchmarkWindow.Reset();
                    }
                }

                sparsePhysics.UpdateDebugLines(sparseRegistry);
                renderSystem.Tick(sparseRegistry, dt, sparsePhysics.GetDebugLines());

                ++frameIndex;
            }

            if (benchmarkWindow.sampleCount > 0) {
                PrintBenchmarkWindow(benchmarkWindow);
            }

            std::cout << "--- EXITING PHYSICS BENCHMARK ---" << std::endl;
            vkDeviceWaitIdle(vulkanDevice.GetNativeDevice());
            renderSystem.Shutdown();
        }

        vkDestroySurfaceKHR(instance, surface, nullptr);
    } catch (const std::exception& e) {
        std::cerr << "\n[FATAL ERROR] " << e.what() << std::endl;
    }

    if (window) {
        glfwDestroyWindow(window);
        glfwTerminate();
    }

    return 0;
}

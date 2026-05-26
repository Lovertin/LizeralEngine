# LizeralEngine

LizeralEngine 是一个个人开发的轻量级现代 3D 引擎。当前工程已经具备可运行编辑器、ECS、Vulkan 渲染、Bullet 物理、资源导入与静态反射生成链路。

这个 README 作为项目根目录的长期入口文档，记录当前架构、构建方式、模块边界和后续演进重点。

## 1. 项目定位

引擎核心方向很明确：

- 以 C++17 为主语言。
- 以 Vulkan 为图形后端。
- 以 Qt6 构建编辑器界面。
- 以纯数据驱动 ECS 组织实体、组件和系统。
- 以 Mesh Shader、Bindless、BDA、硬件光追、SVGF/TAA 等现代 GPU 路线为渲染探索方向。
- 以 Bullet Physics 提供刚体物理。
- 以 Assimp + meshoptimizer 处理模型导入和 meshlet 烘焙。
- 以 Clang/libclang + mustache 模板生成反射和序列化代码。

## 2. 顶层目录

当前主要结构如下：

- `README.md`：引擎理念、ECS 架构、物理系统、相机系统说明。
- `build.ps1` / `build.sh`：构建入口。Windows 下 `build.ps1` 是当前更完整的流程。
- `.tools`：本地构建工具链，包含 CMake、Ninja、Qt、MinGW、VulkanSDK 等路径预设。
- `asset`：测试和演示资源，包括 glb/fbx/hdr、Sponza、车辆、室内场景、天空盒、纹理等。
- `config/reflection_config.json`：反射生成配置。
- `engine/CMakeLists.txt`：主 CMake 工程，定义编辑器和物理沙盒目标。
- `engine/source/runtime`：运行时核心。
- `engine/source/editor`：Qt 编辑器。
- `engine/source/shader`：GLSL/HLSL 风格 shader 源文件和已编译 SPIR-V。
- `engine/source/meta_parser`：静态反射代码生成器。
- `engine/source/_generated`：已生成的 reflection/serializer 代码。
- `engine/bin`：已构建产物和运行时 DLL。

## 3. 构建与运行

Windows 构建脚本是 `D:\LizeralEngine\build.ps1`。它会：

1. 使用 `.tools` 下的 CMake、Ninja、Qt 6.10.0 MinGW、MinGW 13.1、VulkanSDK。
2. 配置 `engine/build-fast`。
3. 构建 Debug 版本。
4. 检查 `engine/bin/LizeralEditor.exe` 和 `engine/bin/LizeralPhysicsSandbox.exe`。
5. 对编辑器执行 `windeployqt`。
6. 默认运行 Editor，也可运行 Sandbox。

常用命令：

```powershell
.\build.ps1 -RunTarget Editor
.\build.ps1 -RunTarget Sandbox
.\build.ps1 -SkipRun
.\build.ps1 -Clean -SkipRun
```

当前 CMake 目标：

- `LizeralEditor`：主编辑器，可视化编辑、视口渲染、场景操作。
- `LizeralPhysicsSandbox`：物理沙盒/测试入口。
- `LizeralRenderingTest`：CMake 中保留了旧目标注释，目前未启用。

## 4. 运行时核心

### 4.1 ECS

路径：`engine/source/runtime/function/ecs`

当前 ECS 有两套实现：

- `Registry`：基础 sparse set ECS。
- `HybridRegistry`：混合 sparse set + archetype/chunk 的实验/优化版本。

基础 `Registry` 的核心机制：

- `Entity` 是轻量整数 ID。
- 每种组件类型对应一个 `Pool<T>`。
- `Pool<T>` 使用 sparse set：
  - `m_components` 连续存储组件数据。
  - `m_dense_entities` 保存 dense index 到 entity 的映射。
  - `m_sparse_indices` 保存 entity 到 dense index 的映射。
  - 删除组件时使用 swap-pop，维持紧凑数组。
- `view<Components...>()` 选择最小组件池作为驱动集合，再过滤拥有全部组件的实体。

`HybridRegistry` 在此基础上增加：

- archetype signature。
- chunk column storage。
- chunk capacity，默认 64。
- archetype query cache。
- 针对标记为 archetype-managed 的组件走 chunk 查询，否则回退 sparse pool。

这说明引擎目前已经从“能用 ECS”进入到“比较不同 ECS 数据布局性能”的阶段。

### 4.2 当前组件集合

路径：`engine/source/runtime/function/ecs/components`

当前主要组件包括：

- `TransformComponent`
- `NameComponent`
- `CameraComponent`
- `CameraControlComponent`
- `ColliderComponent`
- `RigidBodyComponent`
- `DirectionalLightComponent`
- `PointLightComponent`
- `VulkanModelComponent`
- `EditorOnlyComponent`

这套组件覆盖了第一版编辑器和运行时所需的基本对象表达：空间变换、命名、模型、相机、灯光、碰撞、刚体、编辑器专用标记。

### 4.3 物理系统

路径：`engine/source/runtime/function/physics`

物理系统基于 Bullet Physics，核心类是 `PhysicsSystem`。当前能力包括：

- 初始化和关闭物理世界。
- 根据 ECS 组件创建/销毁 Bullet 刚体。
- 每帧 step simulation。
- 将 Bullet 模拟结果写回 `TransformComponent`。
- 支持 `Registry` 与 `HybridRegistry` 两套 ECS。
- 支持 chunked 执行模式。
- 提供 profiling 数据：
  - entity count
  - body create count
  - dirty sync count
  - transform write-back count
  - pre-step/simulation/write-back/total 时间
- 支持 raycast。
- 支持 debug draw line 数据输出。
- 支持 sparse set 与 hybrid archetype 访问性能对比。

整体设计是典型的数据处理系统：物理系统并不理解“物体语义”，只处理带有 `TransformComponent`、`RigidBodyComponent`、`ColliderComponent` 的实体。

### 4.4 相机系统

路径：

- `engine/source/runtime/function/render/CameraSystem`
- `engine/source/runtime/function/render/CameraControlSystem`

相机被拆成：

- `TransformComponent`：空间位置和旋转。
- `CameraComponent`：光学参数和矩阵。
- `CameraControlComponent`：移动速度、鼠标灵敏度等控制参数。

这种拆分让“可控制性”不是 Camera 类的固定能力，而是一个可组合组件。理论上任何具有 Transform 的实体都可以挂上控制组件进行编辑器视口式移动。

## 5. Vulkan 渲染系统

路径：`engine/source/runtime/function/render`

当前渲染模块分层如下：

- `RenderingSystem`：较高层渲染系统入口。
- `VulkanRenderingSystem`：当前主 Vulkan 渲染管线组织者。
- `VulkanRenderer`：渲染器封装。
- `rhi/vulkan`：Vulkan RHI 封装。
- `MeshletBuilder`：meshlet 模型构建相关。

Vulkan RHI 已有封装：

- `VulkanContext`
- `VulkanDevice`
- `VulkanSwapchain`
- `VulkanCommandPool`
- `VulkanCommandBuffer`
- `VulkanBuffer`
- `VulkanTexture`
- `VulkanPipelineBuilder`
- `VulkanDescriptorBuilder`
- `VulkanComputePipeline`
- `VulkanBLAS`
- `VulkanTLAS`

`VulkanRenderingSystem` 当前维护的关键资源：

- 多个 G-Buffer attachment：
  - albedo/metallic
  - normal/roughness
  - depth
  - velocity
  - direct light
  - noisy GI
  - denoised GI
  - GI history
  - moments history
  - final history
- global frame data：
  - viewProj
  - invViewProj
  - prevViewProj
  - camera position
  - light direction/color/intensity
  - frame index
  - jitter
- bindless/global texture list。
- BDA frame/resource buffers。
- per-instance data 和 previous model matrix。
- point light GPU buffer。
- TLAS/BLAS 加速结构。
- 材质 override buffer cache。

当前渲染 preset：

- `Stable`：GI off，hard shadow。
- `SSGI`：screen-space GI，hard shadow。
- `RTGI`：ray-traced GI，soft shadow。

Shader 文件显示当前已经有：

- mesh/task shader：`car.task`、`car.mesh`、`triangle.mesh`
- deferred lighting：`lighting.vert`、`lighting.frag`、`lighting_uber.frag`
- SVGF temporal / a-trous：`svgf_temporal.frag`、`svgf_atrous.frag`
- TAA：`taa.frag`
- transparent pass：`transparent.frag`
- blit：`blit.frag`
- editor grid/axis：`editor_grid.*`、`editor_axis.*`
- debug line：`debug_line.*`
- shared shader library：`common/SceneShared.glsl`、`common/Math.glsl`、`common/BRDF.glsl`、`ssrm_lib.glsl`

## 6. 资源系统

路径：`engine/source/runtime/resource`

当前资源链路大致是：

1. `AssimpModelImporter` 从源模型导入。
2. 生成 `ImportedModelData`：
   - nodes
   - meshes
   - materials
   - material assets
   - encoded textures
3. `MeshletModelCooker` 根据 `MeshletCookOptions` 烘焙为 `RuntimeModelData`。
4. Runtime model 数据包含：
   - mesh assets
   - meshlet vertices
   - micro indices
   - meshlet descriptors
   - bounds
   - textures
   - materials
   - material assets
5. 渲染系统将 runtime model 转成 Vulkan buffer、texture、BLAS 等 GPU 资源。

资源系统已经不是简单“加载一个模型画出来”，而是围绕 meshlet 渲染和 GPU-driven 数据组织在设计。

## 7. 编辑器

路径：`engine/source/editor`

主入口：

- `LizeralEditor.cpp`
- `window/LizeralEditorWindow.*`

`LizeralEditorWindow` 当前持有：

- 全局 `Registry`
- `PhysicsScene`
- `PhysicsSystem`
- `VulkanRenderingSystem`
- `CameraSystem`
- `CameraControlSystem`
- `EngineViewportWidget`
- Scene Outliner
- Inspector
- Control Panel
- engine tick timer
- play mode snapshot

编辑器子系统：

- `viewport/EngineViewportWidget`：Qt 视口承载。
- `panels/SceneOutlinerPanel`：场景实体层级/列表。
- `panels/InspectorPanel`：实体组件查看与编辑。
- `panels/EditorControlPanel`：Play/Stop/Save/Load 工具栏。
- `selection/EditorSelection`：选择状态。
- `context/EditorContext`：编辑器上下文，提供按组件名获取组件等桥接能力。
- `factory/ComponentUIFactory`：组件 UI 生成。
- `factory/DataTypeUIFactory`：基础数据类型 UI 生成。
- `command`：命令系统，含 `ICommand`、`EditPropertyCommand`、`CommandManager`。
- `overlay`：编辑器视口 overlay 数据。
- `style/UE5_Dark.qss`：暗色编辑器风格。

Inspector 当前支持组件注册和 Add Component 菜单，这说明编辑器已经能通过组件组合来驱动实体编辑。

Control Panel 已暴露：

- 进入 Play Mode
- 回到 Edit Mode
- Save Scene
- Load Scene

## 8. 反射与序列化

路径：

- `engine/source/meta_parser`
- `engine/source/runtime/core/meta`
- `engine/source/_generated`
- `engine/template`
- `config/reflection_config.json`

当前反射系统由两部分构成：

- 运行时反射/序列化基础设施：`runtime/core/meta`
- 离线代码生成器：`meta_parser`

`meta_parser` 使用：

- libclang/clang-c 接口读取 C++ AST。
- `__attribute__((annotate(...)))` 风格标注辅助提取类型信息。
- mustache 模板生成 reflection/serializer 文件。

已生成内容包括：

- math 类型：`vector2/3/4`、`quaternion`、`matrix4`、`transform`、`axis_aligned`、`color`
- ECS 组件：`TransformComponent`、`CameraComponent`、`CameraControlComponent`、`ColliderComponent`、`RigidBodyComponent`、`NameComponent`、`ModelComponent`、`VulkanModelComponent`、`DirectionalLightComponent`、`PointLightComponent`
- 资源/材质：`Mesh`、`Material`、`PBRMaterial`、`Texture2D`、`TextureCube`
- 聚合入口：`all_reflection.h`、`all_serializer.h`、`all_serializer.ipp`

当前 `reflection_config.json` 只列出了一部分 math/color 文件，但 `_generated` 目录已经包含更多历史生成结果。后续需要确认生成配置是否已经落后于实际组件集合。

## 9. 第三方依赖

当前主要依赖：

- Qt6：编辑器 UI。
- Vulkan SDK：图形后端。
- GLFW/glad：窗口/GL 辅助历史或测试链路。
- Bullet：物理。
- Assimp：模型导入。
- meshoptimizer：meshlet/网格优化。
- VMA：Vulkan memory allocation。
- stb_image：图片加载。
- json11：JSON 读写。
- libclang/LLVM C API：静态反射解析器。
- mustache.hpp：代码生成模板。

CMake 中针对 Assimp 做了裁剪：

- 关闭 tests/tools。
- 关闭默认全部 importer/exporter。
- 启用 glTF、OBJ、FBX importer。
- 静态构建。

## 10. 当前能力边界

当前已经具备：

- 可构建的 Windows 工具链流程。
- 可运行的 Qt 编辑器目标。
- 可运行的物理沙盒目标。
- ECS 基础架构和 hybrid ECS 实验架构。
- 常见组件集合。
- Bullet 刚体物理同步。
- Vulkan deferred/meshlet/ray tracing 方向的渲染框架。
- G-Buffer、temporal history、SVGF、TAA、debug/editor overlay 等渲染资源布局。
- Assimp 模型导入。
- meshoptimizer meshlet 烘焙。
- 模型、材质、纹理到 GPU 资源的缓存与绑定。
- 编辑器 outliner/inspector/control panel/viewport 基础。
- 反射和序列化生成体系。
- 场景保存/加载相关入口。

尚需进一步确认或补强的地方：

- `reflection_config.json` 与实际 generated 文件范围可能不一致。
- `engine/bin`、`engine/source/_generated`、shader `.spv` 当前都在仓库/目录中，后续需要明确哪些是源资产、哪些是构建产物。
- README 文本存在编码显示问题，疑似 UTF-8 被错误解码，需要统一文档编码。
- CMake 中路径仍硬编码 `C:/Qt/6.10.0/mingw_64`，而 `build.ps1` 使用 `.tools`，两者应继续统一。
- `LizeralRenderingTest` 目标处于注释状态，若后续仍需要，应恢复或移除。
- `DirectionLightComponent` / `DirectionalLightComponent` 命名在部分编辑器上下文里疑似有不一致，需要核对。
- Vendor 体积很大，尤其 Assimp test/model 内容，后续仓库瘦身或 submodule/vendor 策略值得整理。
- 暂未看到独立自动化测试体系，后续至少可补 ECS、资源 cooker、序列化、物理同步的单元测试。

## 11. 后续文档规划

建议后续继续补充：

- `docs/Architecture`：模块架构、数据流、渲染管线图。
- `docs/Stage_Summary`：阶段总结。
- `docs/Exec_Plans`：每次较大改动的执行计划。
- `docs/Research`：Vulkan、mesh shader、RTGI、SVGF、ECS 性能实验记录。
- `docs/Issues`：当前风险、bug、技术债和待验证清单。

## 12. 当前整体判断

LizeralEngine 当前的核心价值不在“功能铺满”，而在它已经搭起了一个方向一致的现代小型引擎骨架：

- 数据层用 ECS 维持清晰边界。
- 资源层面向 meshlet 和 GPU 数据布局。
- 渲染层直接拥抱 Vulkan、BDA、mesh shader、RT/TAA/SVGF。
- 编辑器通过 Qt 和反射系统把组件暴露出来。
- 物理系统作为 ECS processor 独立同步。

这已经是一个可以继续演进的技术底座。接下来最值得优先做的是：统一生成链路、收敛构建产物策略、补关键测试、把渲染/资源/编辑器的数据流文档化，然后再围绕一个真实 demo 场景持续打磨稳定性。

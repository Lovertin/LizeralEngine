#pragma once

#include "runtime/function/physics/PhysicsSystem.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace Lizeral {

    enum class EditorOverlayLineLayer : uint8_t {
        Physics = 0,
        Selection,
        CameraFrustum,
        Gizmo,
        Custom,
        Count
    };

    constexpr size_t EditorOverlayLineLayerCount = static_cast<size_t>(EditorOverlayLineLayer::Count);

    struct EditorOverlayGridSettings {
        bool enabled { true };
        float planeHeight { 0.0f };
        float minorSpacing { 1.0f };
        float majorSpacing { 10.0f };
        float fadeDistance { 160.0f };
        float minorOpacity { 0.20f };
        float majorOpacity { 0.45f };
        float axisOpacity { 0.95f };
    };

    struct EditorOverlayAxisSettings {
        bool showAxes { true };
        float axisHalfExtent { 100.0f };
    };

    struct EditorViewportOverlayData {
        bool enabled { true };
        EditorOverlayGridSettings grid {};
        EditorOverlayAxisSettings axes {};
        std::vector<DebugLineVertex> worldLines;

        bool HasVisibleContent() const {
            return enabled && (grid.enabled || axes.showAxes || !worldLines.empty());
        }
    };

} // namespace Lizeral

#pragma once

#include "editor/overlay/EditorViewportOverlayTypes.h"

#include <array>

namespace Lizeral {

    class EditorViewportOverlay {
    public:
        EditorViewportOverlay() = default;

        void SetEnabled(bool enabled);

        EditorOverlayGridSettings& GridSettings() { return m_frameData.grid; }
        const EditorOverlayGridSettings& GridSettings() const { return m_frameData.grid; }

        void SetGridEnabled(bool enabled);

        void SetPhysicsDebugLines(const std::vector<DebugLineVertex>& lines);
        void SetLineLayer(EditorOverlayLineLayer layer, const std::vector<DebugLineVertex>& lines);
        std::vector<DebugLineVertex>& AccessLineLayer(EditorOverlayLineLayer layer);
        const std::vector<DebugLineVertex>& GetLineLayer(EditorOverlayLineLayer layer) const;

        void ClearLineLayer(EditorOverlayLineLayer layer);
        void ClearRuntimeToolLayers();

        const EditorViewportOverlayData& BuildFrameData();

    private:
        static size_t ToIndex(EditorOverlayLineLayer layer);
        void RebuildMergedLines();

        bool m_dirty { true };
        EditorViewportOverlayData m_frameData {};
        std::array<std::vector<DebugLineVertex>, EditorOverlayLineLayerCount> m_lineLayers;
    };

} // namespace Lizeral

#include "EditorViewportOverlay.h"

#include <algorithm>

namespace Lizeral {

    size_t EditorViewportOverlay::ToIndex(EditorOverlayLineLayer layer) {
        return static_cast<size_t>(layer);
    }

    void EditorViewportOverlay::SetEnabled(bool enabled) {
        if (m_frameData.enabled == enabled) {
            return;
        }

        m_frameData.enabled = enabled;
        m_dirty = true;
    }

    void EditorViewportOverlay::SetGridEnabled(bool enabled) {
        if (m_frameData.grid.enabled == enabled) {
            return;
        }

        m_frameData.grid.enabled = enabled;
        m_dirty = true;
    }

    void EditorViewportOverlay::SetPhysicsDebugLines(const std::vector<DebugLineVertex>& lines) {
        SetLineLayer(EditorOverlayLineLayer::Physics, lines);
    }

    void EditorViewportOverlay::SetLineLayer(EditorOverlayLineLayer layer, const std::vector<DebugLineVertex>& lines) {
        auto& target = m_lineLayers[ToIndex(layer)];
        target = lines;
        m_dirty = true;
    }

    std::vector<DebugLineVertex>& EditorViewportOverlay::AccessLineLayer(EditorOverlayLineLayer layer) {
        m_dirty = true;
        return m_lineLayers[ToIndex(layer)];
    }

    const std::vector<DebugLineVertex>& EditorViewportOverlay::GetLineLayer(EditorOverlayLineLayer layer) const {
        return m_lineLayers[ToIndex(layer)];
    }

    void EditorViewportOverlay::ClearLineLayer(EditorOverlayLineLayer layer) {
        auto& target = m_lineLayers[ToIndex(layer)];
        if (target.empty()) {
            return;
        }

        target.clear();
        m_dirty = true;
    }

    void EditorViewportOverlay::ClearRuntimeToolLayers() {
        for (size_t i = ToIndex(EditorOverlayLineLayer::Selection); i < m_lineLayers.size(); ++i) {
            if (!m_lineLayers[i].empty()) {
                m_lineLayers[i].clear();
                m_dirty = true;
            }
        }
    }

    const EditorViewportOverlayData& EditorViewportOverlay::BuildFrameData() {
        if (m_dirty) {
            RebuildMergedLines();
        }

        return m_frameData;
    }

    void EditorViewportOverlay::RebuildMergedLines() {
        size_t totalVertexCount = 0;
        for (const auto& layer : m_lineLayers) {
            totalVertexCount += layer.size();
        }

        m_frameData.worldLines.clear();
        m_frameData.worldLines.reserve(totalVertexCount);

        for (const auto& layer : m_lineLayers) {
            m_frameData.worldLines.insert(m_frameData.worldLines.end(), layer.begin(), layer.end());
        }

        m_dirty = false;
    }

} // namespace Lizeral

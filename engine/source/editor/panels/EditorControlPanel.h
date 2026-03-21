#pragma once
#include <QToolBar>
#include <QPushButton>

namespace Lizeral {

    class EditorControlPanel : public QToolBar {
        Q_OBJECT
    public:
        explicit EditorControlPanel(QWidget* parent = nullptr);

    signals:
        void OnPlayModeEntered();
        void OnEditModeEntered();

        void OnSaveScene();
        void OnLoadScene();

    private:
        QPushButton* m_PlayButton { nullptr };
        QPushButton* m_StopButton { nullptr };

        QPushButton* m_SaveButton { nullptr };
        QPushButton* m_LoadButton { nullptr };
    };

} // namespace Lizeral
#pragma once
#include <QToolBar>
#include <QPushButton>

namespace Lizeral {

    class EditorControlPanel : public QToolBar {
        Q_OBJECT
    public:
        explicit EditorControlPanel(QWidget* parent = nullptr);

    signals:
        // 向外广播状态切换事件，让 MainWindow 去处理具体的场景序列化和面板刷新
        void OnPlayModeEntered();
        void OnEditModeEntered();

    private:
        QPushButton* m_PlayButton { nullptr };
        QPushButton* m_StopButton { nullptr };
    };

} // namespace Lizeral
#include "EditorControlPanel.h"
#include "editor/context/EditorContext.h"
#include <iostream>

namespace Lizeral {

    EditorControlPanel::EditorControlPanel(QWidget* parent) : QToolBar(parent) {
        setMovable(false); // 固定在顶部，不让用户拖拽成悬浮窗
        setStyleSheet("QToolBar { background-color: #333333; border-bottom: 1px solid #222; padding: 4px; }");

        // 1. 创建按钮
        m_PlayButton = new QPushButton("▶ Play", this);
        m_StopButton = new QPushButton("■ Stop", this);

        // 设置样式
        m_PlayButton->setStyleSheet("QPushButton { background-color: #4CAF50; color: white; font-weight: bold; padding: 6px 20px; border-radius: 3px; margin: 2px; } QPushButton:hover { background-color: #45a049; }");
        m_StopButton->setStyleSheet("QPushButton { background-color: #c9302c; color: white; font-weight: bold; padding: 6px 20px; border-radius: 3px; margin: 2px; } QPushButton:hover { background-color: #d9534f; }");
        
        // 初始状态：只能点 Play
        m_StopButton->setEnabled(false);

        // 2. 添加到工具栏 (你可以加一些 Spacer 让按钮居中，这里先简单按顺序添加)
        addWidget(m_PlayButton);
        addWidget(m_StopButton);

        // 3. 绑定点击逻辑
        connect(m_PlayButton, &QPushButton::clicked, this, [this]() {
            m_PlayButton->setEnabled(false);
            m_StopButton->setEnabled(true);
            
            // 进入 Play 模式前，清空撤销栈（防止 Play 运行到一半按 Ctrl+Z 搞乱数据）
            EditorContext::Get().GetCommandManager()->Clear();
            
            // 切换全局状态机
            EditorContext::Get().setMode(EditorMode::Play);
            
            // 抛出信号，让 MainWindow 去做场景快照
            emit OnPlayModeEntered();
        });

        connect(m_StopButton, &QPushButton::clicked, this, [this]() {
            m_PlayButton->setEnabled(true);
            m_StopButton->setEnabled(false);
            
            // 切回 Edit 模式
            EditorContext::Get().setMode(EditorMode::Edit);
            
            // 抛出信号，让 MainWindow 去恢复场景快照并刷新 UI
            emit OnEditModeEntered();
        });
    }

} // namespace Lizeral
#include "EditorControlPanel.h"
#include "editor/context/EditorContext.h"
#include <iostream>
#include <QWidget> 

namespace Lizeral {

    EditorControlPanel::EditorControlPanel(QWidget* parent) : QToolBar(parent) {
        setMovable(false); 
        setStyleSheet("QToolBar { background-color: #333333; border-bottom: 1px solid #222; padding: 4px; }");

        m_SaveButton = new QPushButton("💾 Save", this);
        m_LoadButton = new QPushButton("📂 Load", this);

        QString normalStyle = "QPushButton { background-color: #2b2b2b; color: white; padding: 6px 15px; border-radius: 3px; margin: 2px; border: 1px solid #444; } QPushButton:hover { background-color: #444; }";
        m_SaveButton->setStyleSheet(normalStyle);
        m_LoadButton->setStyleSheet(normalStyle);

        addWidget(m_SaveButton);
        addWidget(m_LoadButton);

        QWidget* leftSpacer = new QWidget(this);
        leftSpacer->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Preferred);
        addWidget(leftSpacer);

        //Core: Single Play/Stop toggle button (centered)
        m_PlayButton = new QPushButton("▶ Play", this); 

        QString playStyle = "QPushButton { background-color: #4CAF50; color: white; font-weight: bold; padding: 6px 20px; border-radius: 3px; margin: 2px; } QPushButton:hover { background-color: #45a049; }";
        QString stopStyle = "QPushButton { background-color: #c9302c; color: white; font-weight: bold; padding: 6px 20px; border-radius: 3px; margin: 2px; } QPushButton:hover { background-color: #d9534f; }";
        
        m_PlayButton->setStyleSheet(playStyle);
        addWidget(m_PlayButton);

        QWidget* rightSpacer = new QWidget(this);
        rightSpacer->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Preferred);
        addWidget(rightSpacer);

        // bind click logic
        connect(m_SaveButton, &QPushButton::clicked, this, [this]() { emit OnSaveScene(); });
        connect(m_LoadButton, &QPushButton::clicked, this, [this]() { emit OnLoadScene(); });

        connect(m_PlayButton, &QPushButton::clicked, this, [this, playStyle, stopStyle]() mutable {
            
            if (m_PlayButton->text() == "▶ Play") {
                m_PlayButton->setText("■ Stop");
                m_PlayButton->setStyleSheet(stopStyle);
                
                // swap
                EditorContext::Get().GetCommandManager()->Clear();
                EditorContext::Get().setMode(EditorMode::Play);
                emit OnPlayModeEntered();
            } 
            else {
                m_PlayButton->setText("▶ Play");
                m_PlayButton->setStyleSheet(playStyle);
       
                //swap
                EditorContext::Get().setMode(EditorMode::Edit);
                emit OnEditModeEntered();
            }
        });
    }

} // namespace Lizeral
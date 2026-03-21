#include "CommandManager.h"

namespace Lizeral {

    void CommandManager::ExecuteCommand(std::unique_ptr<ICommand> command) {
        if (!command) return;

        if (!m_UndoStack.empty()) {
            if (m_UndoStack.back()->MergeWith(command.get())) {
                // merge success
                m_UndoStack.back()->Execute();
                return;
            }
        }

        // push stack
        command->Execute();
        m_UndoStack.push_back(std::move(command));

        // clear Redo stack
        m_RedoStack.clear();

        if (m_UndoStack.size() > m_MaxHistory) {
            m_UndoStack.pop_front(); 
        }
    }

    void CommandManager::Undo() {
        if (m_UndoStack.empty()) return;

        auto cmd = std::move(m_UndoStack.back());
        m_UndoStack.pop_back();

        cmd->Undo();

        m_RedoStack.push_back(std::move(cmd));
    }

    void CommandManager::Redo() {
        if (m_RedoStack.empty()) return;

        auto cmd = std::move(m_RedoStack.back());
        m_RedoStack.pop_back();

        cmd->Execute();

        m_UndoStack.push_back(std::move(cmd));
    }

    void CommandManager::Clear() {
        m_UndoStack.clear();
        m_RedoStack.clear();
    }

} // namespace Lizeral
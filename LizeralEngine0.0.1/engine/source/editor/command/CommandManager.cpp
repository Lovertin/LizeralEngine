#include "CommandManager.h"

namespace Lizeral {

    void CommandManager::ExecuteCommand(std::unique_ptr<ICommand> command) {
        if (!command) return;

        // 1. 尝试与栈顶的指令合并 (比如连续按了几十次 SpinBox 的上箭头)
        if (!m_UndoStack.empty()) {
            if (m_UndoStack.back()->MergeWith(command.get())) {
                // 合并成功，直接执行栈顶指令（此时它的新值已经被更新）
                m_UndoStack.back()->Execute();
                return;
            }
        }

        // 2. 无法合并，作为新指令执行并压栈
        command->Execute();
        m_UndoStack.push_back(std::move(command));

        // 3. 时间线产生分支，必须清空 Redo 栈
        m_RedoStack.clear();

        // 4. 堆栈容量控制
        if (m_UndoStack.size() > m_MaxHistory) {
            m_UndoStack.pop_front(); 
        }
    }

    void CommandManager::Undo() {
        if (m_UndoStack.empty()) return;

        auto cmd = std::move(m_UndoStack.back());
        m_UndoStack.pop_back();

        cmd->Undo(); // 内部会安全地重新获取组件指针并写入旧值

        m_RedoStack.push_back(std::move(cmd));
    }

    void CommandManager::Redo() {
        if (m_RedoStack.empty()) return;

        auto cmd = std::move(m_RedoStack.back());
        m_RedoStack.pop_back();

        cmd->Execute(); // 写入新值

        m_UndoStack.push_back(std::move(cmd));
    }

    void CommandManager::Clear() {
        m_UndoStack.clear();
        m_RedoStack.clear();
    }

} // namespace Lizeral
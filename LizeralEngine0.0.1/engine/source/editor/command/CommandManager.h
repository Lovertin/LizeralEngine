#pragma once
#include "ICommand.h"
#include <deque>
#include <memory>

namespace Lizeral {

    class CommandManager {
    public:
        // 限制最大历史记录为 100 步，防止内存溢出
        CommandManager(size_t maxHistory = 100) : m_MaxHistory(maxHistory) {}

        void ExecuteCommand(std::unique_ptr<ICommand> command);
        void Undo();
        void Redo();
        void Clear();

        size_t GetUndoSize() const { return m_UndoStack.size(); }
        size_t GetRedoSize() const { return m_RedoStack.size(); }

    private:
        size_t m_MaxHistory;
        std::deque<std::unique_ptr<ICommand>> m_UndoStack; 
        std::deque<std::unique_ptr<ICommand>> m_RedoStack;
    };

} // namespace Lizeral
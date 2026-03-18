#pragma once
#include "ICommand.h"
#include <deque>
#include <memory>

namespace Lizeral {

    class CommandManager {
    public:
        // 限制最大历史记录为 50 步，防止内存溢出
        CommandManager(size_t maxHistory = 50) : m_MaxHistory(maxHistory) {}

        void ExecuteCommand(std::unique_ptr<ICommand> command);
        void Undo();
        void Redo();
        void Clear();

    private:
        size_t m_MaxHistory;
        std::deque<std::unique_ptr<ICommand>> m_UndoStack; 
        std::deque<std::unique_ptr<ICommand>> m_RedoStack;
    };

} // namespace Lizeral
#pragma once

namespace Lizeral {

    class ICommand {
    public:
        virtual ~ICommand() = default;

        // 正向执行（或重做时调用）
        virtual void Execute() = 0;

        // 撤销时调用
        virtual void Undo() = 0;

        // 【核心机制】：尝试与上一个指令合并
        // 用于处理高频的连续操作（比如连续点击按钮数十次）
        virtual bool MergeWith(const ICommand* other) { 
            return false; // 默认不允许合并
        }
    };

} // namespace Lizeral
#pragma once

namespace Lizeral {

    class ICommand {
    public:
        virtual ~ICommand() = default;

        virtual void Execute() = 0;

        virtual void Undo() = 0;

        virtual bool MergeWith(const ICommand* other) { 
            return false;
        }
    };

} // namespace Lizeral
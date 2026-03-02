#include <assert.h>

#include "runtime/core/meta/json.h"
#include "runtime/core/meta/serializer/serializer.h"
#include "runtime/core/meta/reflection/reflection.h"
#include "runtime/core/meta/reflection/reflection_register.h"

#include "_generated/reflection/all_reflection.h"
#include "_generated/serializer/all_serializer.h"
#include "_generated/serializer/all_serializer.ipp"

namespace Lizeral
{
    namespace Reflection
    {
        void TypeMetaRegister::Unregister() { TypeMetaRegisterinterface::unregisterAll(); }
    } // namespace Reflection
} // namespace Piccolo

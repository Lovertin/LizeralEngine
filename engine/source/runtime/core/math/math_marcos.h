#pragma once

#define LIZERAL_MIN(x, y) (((x) < (y)) ? (x) : (y))
#define LIZERAL_MAX(x, y) (((x) > (y)) ? (x) : (y))
#define LIZERAL_PIN(a, min_value, max_value) LIZERAL_MIN(max_value, LIZERAL_MAX(a, min_value))

#define LIZERAL_VALID_INDEX(idx, range) (((idx) >= 0) && ((idx) < (range)))
#define LIZERAL_PIN_INDEX(idx, range) LIZERAL_PIN(idx, 0, (range)-1)

#define LIZERAL_SIGN(x) ((((x) > 0.0f) ? 1.0f : 0.0f) + (((x) < 0.0f) ? -1.0f : 0.0f))

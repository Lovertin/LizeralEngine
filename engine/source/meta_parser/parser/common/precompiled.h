#pragma once

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <clang-c/Index.h>

namespace fs = std::filesystem;

#include "utils/meta_data_config.h"
#include "utils/meta_utils.h"

#include "mustache.hpp"

namespace Mustache = kainjow::mustache;

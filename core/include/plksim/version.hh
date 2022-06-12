#pragma once

namespace plksim {

#define MAJOR_VERSION 0
#define MINOR_VERSION 1
#define PATCH_VERSION 0

#define STR(x) STR_HELPER(x)
#define STR_HELPER(x) #x
#define VERSION STR(MAJOR_VERSION) "." STR(MINOR_VERSION) "." STR(PATCH_VERSION)

} // namespace plksim
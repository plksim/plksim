cmake_minimum_required(VERSION 3.20)
project("plksim")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

enable_testing()

message("CMAKE_SYSTEM_NAME=${CMAKE_SYSTEM_NAME}")
message("CMAKE_SYSTEM_PROCESSOR=${CMAKE_SYSTEM_PROCESSOR}")
message("CMAKE_SYSTEM_NAME=${CMAKE_SYSTEM_NAME}")


# =============================================================================
# deal.II
# =============================================================================
if (APPLE)
  set(DEAL_II_DIR /Applications/deal.II.app/Contents/Resources/Libraries/lib/cmake/deal.II)
endif()
find_package(deal.II REQUIRED
    HINTS ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
    )
deal_ii_initialize_cached_variables()
message("deal.II_VERSION=${DEAL_II_VERSION}")


add_subdirectory(${PROJECT_SOURCE_DIR}/core)

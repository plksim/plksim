# =============================================================================
# Library
# =============================================================================
file(GLOB_RECURSE PUBLIC_HEADERS "include/**/*.hh")
file(GLOB_RECURSE PRIVATE_HEADERS "src/**/*.hh")
file(GLOB_RECURSE SRCS "src/**/*.cc")

set(LIBS
    )

add_library(plksim ${PUBLIC_HEADERS} ${PRIVATE_HEADERS} ${SRCS})
deal_ii_setup_target(plksim)
target_link_libraries(plksim PRIVATE ${LIBS})

install(TARGETS plksim 
    ARCHIVE DESTINATION lib
    )
install(DIRECTORY
    "include/plksim" DESTINATION include
)


# =============================================================================
# Tests
# =============================================================================
find_package(GTest REQUIRED)
message("GTest_VERSION=${GTest_VERSION}")

file(GLOB_RECURSE TESTS "test/**/*.test.cc")

add_executable(plksim_test ${TESTS})
target_link_libraries(plksim_test PRIVATE plksim ${LIBS})
target_link_libraries(plksim_test PRIVATE GTest::gtest GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(plksim_test)

# =============================================================================
# Library
# =============================================================================
set(HEADERS
	include/plksim/version.hh
	src/utils/timer.hh
	)

set(SRCS
	src/utils/timer.cc
	)

set(LIBS
	)

add_library(plksim ${HEADERS} ${SRCS})
deal_ii_setup_target(plksim)
target_link_libraries(plksim PRIVATE ${LIBS})


# =============================================================================
# Tests
# =============================================================================
find_package(GTest REQUIRED)
message("GTest_VERSION=${GTest_VERSION}")

set(TESTS
	test/utils/timer.test.cc
	test/version.test.cc
	)

add_executable(plksim_test ${TESTS})
target_link_libraries(plksim_test PRIVATE plksim ${LIBS})
target_link_libraries(plksim_test PRIVATE GTest::gtest GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(plksim_test)

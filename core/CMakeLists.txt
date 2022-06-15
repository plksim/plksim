# =============================================================================
# Sources
# =============================================================================
file(GLOB_RECURSE PUBLIC_HEADERS "include/*.hh")
file(GLOB_RECURSE PRIVATE_HEADERS "src/*.hh")
file(GLOB_RECURSE SRCS "src/*.cc")


# =============================================================================
# 3rd libraries
# =============================================================================
if (APPLE)
    set(DEAL_II_DIR /Applications/deal.II.app/Contents/Resources/Libraries/lib/cmake/deal.II)
endif()
find_package(deal.II REQUIRED
    HINTS ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
    )
deal_ii_initialize_cached_variables()
message("deal.II_VERSION=${DEAL_II_VERSION}")

set(LIBS
    )


# =============================================================================
# define library
# =============================================================================
add_library(${PROJECT_NAME} ${PUBLIC_HEADERS} ${PRIVATE_HEADERS} ${SRCS})
deal_ii_setup_target(${PROJECT_NAME})
target_link_libraries(${PROJECT_NAME} ${LIBS})

include(CMakePackageConfigHelpers)
set(LIB_VERSION ${PROJECT_VERSION})
set(LIB_LIBRARIES ${PROJECT_NAME} ${DEAL_II_LIBRARIES_RELEASE} ${LIBS})
configure_file(include/${PROJECT_NAME}/version.hh.cmake.in
    ${CMAKE_CURRENT_SOURCE_DIR}/include/${PROJECT_NAME}/version.hh)
configure_package_config_file(config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
    INSTALL_DESTINATION lib/cmake/${PROJECT_NAME}
    )
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMinorVersion
    )


# =============================================================================
# define install
# =============================================================================
install(TARGETS ${PROJECT_NAME}
    ARCHIVE DESTINATION lib
    )
install(DIRECTORY include/${PROJECT_NAME}
    DESTINATION include
    )
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
    DESTINATION lib/cmake/${PROJECT_NAME}
    )


# =============================================================================
# Tests
# =============================================================================
find_package(GTest REQUIRED)
message("GTest_VERSION=${GTest_VERSION}")

file(GLOB_RECURSE TESTS "test/*.test.cc")

add_executable(${PROJECT_NAME}_test ${TESTS})
target_link_libraries(${PROJECT_NAME}_test PRIVATE ${PROJECT_NAME} ${LIBS})
target_link_libraries(${PROJECT_NAME}_test PRIVATE GTest::gtest GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(${PROJECT_NAME}_test)

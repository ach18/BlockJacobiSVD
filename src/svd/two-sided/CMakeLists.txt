set(SOURCES
    svd_blocked.cpp
    svd_subprocedure.cpp
    nsvd.cpp
    )

add_library(two-sided-svd-lib STATIC ${SOURCES})
target_include_directories(two-sided-svd-lib PUBLIC .)
target_link_libraries(two-sided-svd-lib utils-lib)

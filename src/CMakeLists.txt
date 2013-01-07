##TODO afegir altres fitxers.
# message(STATUS "C_flags: ${CMAKE_C_FLAGS}\n")
# set(CMAKE_C_FLAGS "-std=gnu99") 
add_library(probability SHARED betabinom.c betanbinom.c beta.c)
add_library(models SHARED cohortModel.c cohortModel-R.c)
message(STATUS "R: ${LIBR_LIBRARIES}\n")

add_executable(lhr main.c)
target_link_libraries(probability ${LIBR_LIBRARIES}) 

target_link_libraries(lhr probability models ${LIBR_LIBRARIES}) 

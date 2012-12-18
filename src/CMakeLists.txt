##TODO afegir altres fitxers.
add_library(probability SHARED betabinom.c betanbinom.c beta.c)
add_library(models SHARED cohortModel.c)
message(STATUS "R: ${LIBR_LIBRARIES}\n")

add_executable(lhr main.c)
target_link_libraries(probability ${LIBR_LIBRARIES}) 

target_link_libraries(lhr probability models ${LIBR_LIBRARIES}) 

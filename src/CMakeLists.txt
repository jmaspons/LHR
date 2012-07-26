
##TODO afegir altres fitxers.
add_library(probability SHARED betabinom.c betanbinom.c)
# add_library(models SHARED cohortModel.c)

# add_executable(lhr main.cpp)
message(STATUS "R: ${LIBR_LIBRARIES}\n")
# target_link_libraries(lhr probability model ${LIBR_LIBRARIES}) 
target_link_libraries(probability ${LIBR_LIBRARIES}) 


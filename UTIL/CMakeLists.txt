set (CMAKE_Fortran_FLAGS "-fixed -extend-source 132")

set(thermo_code SRC/cape_analysis.f SRC/hbuf_lib.f SRC/sat.f SRC/cape.f SRC/cin.f)
add_library(thermo STATIC ${thermo_code})

add_executable(stat2nc SRC/stat2nc.f)
target_link_libraries(stat2nc ${netcdff_lib} ${netcdf_lib} thermo)

add_executable(bin3D2nc SRC/bin3D2nc.f)
target_link_libraries(bin3D2nc ${netcdff_lib} ${netcdf_lib} thermo)

set(executables stat2nc bin3D2nc)
set(libraries thermo)

install(TARGETS ${executables} ${libraries} RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/bin ARCHIVE DESTINATION ${CMAKE_SOURCE_DIR}/lib)

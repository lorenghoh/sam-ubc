cmake_minimum_required(VERSION 3.6)
project(sam_cmake)

# Edit: Link external libraries necessary
set(nc_dir "/usr/local/netcdf/4.6.1_intel18")
set(ncf_dir "/usr/local/netcdf/4.6.1_intel18")
set(mpi_dir "/usr/local/mpi/intel18/openmpi-2.0.2")

enable_language(Fortran)

INCLUDE_DIRECTORIES(${CMAKE_Fortran_MODULE_DIRECTORY} ${nc_dir}/include ${ncf_dir}/include ${mpi_dir}/include ${CMAKE_SOURCE_DIR}/SRC/TIMING)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

# Edit: C and Fortran flags
set(FFLAGS "-O2 -fPIC -r8")
set(CMAKE_Fortran_FLAGS "-free -extend-source 132 ${FFLAGS}")
set(CMAKE_F77_FLAGS "-fixed ${FFLAGS}")

set(CMAKE_Fortran_MODULE_DIRECTORY  ${PROJECT_BINARY_DIR}/include_fortran)

set(LIBRARY_OUTPUT_PATH  ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR/bin})

# Strictly necessary files
FILE(GLOB adv_src "SRC/ADV_MPDATA/*.f90")
list(FILTER adv_src EXCLUDE REGEX ".*/[._#].*.f90")
set(
	first_pass_list
		SRC/domain.f90
		SRC/grid.f90
		SRC/params.f90
		SRC/vars.f90
		SRC/sat.f90
		SRC/tracers.f90
		SRC/hbuffer.f90
		SRC/stat_sw2.f90
		SRC/stat_varscalar.f90
)
add_library(first_pass STATIC ${first_pass_list} ${adv_src})

# Main SAM source files
FILE(GLOB f90_src "SRC/*.f90")
FILE(GLOB f77_src "SRC/*.f")
set_source_files_properties(${f77_src} COMPILE_FLAGS "${CMAKE_F77_FLAGS}")
LIST(REMOVE_ITEM f90_src ${first_pass_list} "SRC/main.f90")
add_library(sam_src STATIC ${f90_src} ${f77_src})

# Dynamic model components
# Microphysics

## SAM 1-moment microphysics scheme
# FILE(GLOB micro_src "SRC/MICRO_SAM1MOM/*.f90")
# list(APPEND micro_src "SRC/precip_fall.f90")

## Morrison's 2-moment microphysics scheme
# FILE(GLOB micro_src "SRC/MICRO_M2005/*.f90")

## Thomson microphysics scheme
FILE(GLOB micro_src "SRC/MICRO_THOM/*.f90")
list(FILTER micro_src EXCLUDE REGEX ".*/[._#].*.f90")

# Radiation
FILE(GLOB rad_src "SRC/RAD_CAM/*.f90")
list(FILTER rad_src EXCLUDE REGEX ".*/[._#].*.f90")

# Entrainment
FILE(GLOB ent_src "SRC/ENTRAINMENT/*.f90")

# SGS TKE
FILE(GLOB sgs_src "SRC/SGS_TKE/*.f90")
list(FILTER sgs_src EXCLUDE REGEX ".*/[._#].*.f90")

# SLM
FILE(GLOB slm_src "SRC/SLM/*.f90")
list(FILTER slm_src EXCLUDE REGEX ".*/[._#].#.f90")

add_library(
	dyn_src STATIC 
		${micro_src}
		${rad_src}
		${ent_src}
		${sgs_src}
		${slm_src}
)

FILE(GLOB sim_f90 "SRC/SIMULATORS/*.f90")
FILE(GLOB sim_f77 "SRC/SIMULATORS/*.f")
set_source_files_properties(${sim_f77} PROPERTIES COMPILE_FLAGS "${CMAKE_F77_FLAGS}")
add_library(sim_src STATIC ${sim_f90} ${sim_f77})
message("Simulators consists of ${sim_f90} ${sim_f77}")

add_dependencies(dyn_src first_pass)
add_dependencies(sam_src first_pass dyn_src sim_src)
add_dependencies(sim_src first_pass dyn_src)

add_executable(SAM_CMAKE "SRC/main.f90")

find_library(mpi_lib mpi ${mpi_dir}/lib NO_DEFAULT_PATH)
find_library(netcdf_lib NAMES netcdf PATHS ${nc_dir}/lib NO_DEFAULT_PATH)
find_library(netcdff_lib NAMES netcdff PATHS ${ncf_dir}/lib NO_DEFAULT_PATH)

target_link_libraries(
	SAM_CMAKE
		sam_src
		sim_src
		dyn_src
		first_pass
		${netcdf_lib}
		${netcdff_lib}
)

add_subdirectory(UTIL)

set(target "${CMAKE_SOURCE_DIR}/RUNDATA")
set(link "${CMAKE_BINARY_DIR}/RUNDATA")
add_custom_target(
	link_rundata ALL
		  COMMAND ${CMAKE_COMMAND} -E create_symlink ${target} ${link}
)

install(TARGETS SAM_CMAKE RUNTIME DESTINATION ${PROJECT_BINARY_DIR}/bin)

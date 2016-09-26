cmake_minimum_required(VERSION 2.6)

project(umat-lib Fortran)

# list of source files in current directory
set (toplevelFiles
		aux_routines.f90
		constants.f90
		derived_types.f90
		interface_definitions.f90
		material_info.f90
	)
	
# collect source files from ./return_mapping
include (return_mapping/CMakeLists.txt)
foreach (FILE ${FILES})
	set (return_mappingFiles ${return_mappingFiles} return_mapping/${FILE})
endforeach (FILE)

# collect source files from ./return_mapping/model_components_all
include (return_mapping/model_components_all/CMakeLists.txt)
foreach (FILE ${FILES})
	set (model_componentsFiles ${model_componentsFiles} return_mapping/model_components_all/${FILE})
endforeach (FILE)

# create static library
add_library(UmatFortran ${toplevelFiles} ${return_mappingFiles} ${model_componentsFiles} )
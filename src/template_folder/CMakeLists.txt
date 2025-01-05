# Discover and install all header files
file(GLOB installable_headers "*.h")

# Get the relative path from the root
file(RELATIVE_PATH relative_path_from_root ${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

# Remove the 'src/' prefix from the relative path if it exists
string(REPLACE "src/" "" relative_path_without_src ${relative_path_from_root})

# Install the header files to the relative path
install(FILES ${installable_headers} DESTINATION include/${project_name}/${relative_path_from_root})

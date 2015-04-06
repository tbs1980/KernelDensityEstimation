find_path(NANOFLANN_INCLUDE_DIR
  NAMES
  nanoflann.hpp
  PATHS
  ${NANOFLANN_ROOT}
  ${INCLUDE_INSTALL_DIR}
  PATH_SUFFIXES
  include
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(nanoflann DEFAULT_MSG NANOFLANN_INCLUDE_DIR)
mark_as_advanced(NANOFLANN_INCLUDE_DIR)

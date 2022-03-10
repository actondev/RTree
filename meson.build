project('RTree', 'cpp')

include_dirs = include_directories('.')

RTree_dep = declare_dependency(
  include_directories: include_dirs
)

if get_option('test')
  subdir('test')
endif

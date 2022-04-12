project('drtree', 'cpp',
       	default_options: [
	  'cpp_std=c++17',
	],
       )

include_dirs = include_directories('.')

drtree_dep = declare_dependency(
  include_directories: include_dirs,
  sources: [
    'aod/rtree_base.cpp',
    ]
)

if get_option('test')
  subdir('test')
endif

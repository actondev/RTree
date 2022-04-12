project('drtree', 'cpp',
       	default_options: [
	  'cpp_std=c++17',
	],
       )

include_dirs = include_directories('.')

rtree = library('rtree',
                sources: [
                  'aod/rtree_base.cpp',
                ],
                include_directories:  [
                  include_dirs,
                  # include_directories('./test') # debugging
                  ]
)

rtree_dep = declare_dependency(
  include_directories: include_dirs,
  link_with: rtree
)


if get_option('test')
  subdir('test')
endif

project('drtree', 'cpp',
       	default_options: [
	  'cpp_std=c++17',
	],
       )

include_dirs = include_directories('.')

rtree = static_library('rtree',
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
  # performance impact when using linked version?
  # link_with: rtree,
  sources: [
    'aod/rtree_base.cpp',
  ],
)


if get_option('test')
  subdir('test')
endif

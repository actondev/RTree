catch2_proj = subproject('catch2')

catch2_dep = catch2_proj.get_variable('catch2_dep')


test_dep = declare_dependency(
  include_directories: include_directories('.')
)

test_src = files( 'main.cpp',
                  'test.cpp'
                )

test_deps = [ declare_dependency(
  include_directories: include_directories('.')  ) ]

run_tests = executable(
  'run-tests',
  test_src,
  dependencies : [
    catch2_dep,
    rtree_dep,
  ],
  cpp_args: [
    '-Wno-narrowing',
    '-Wno-sign-compare'
    # '-DCATCH_CONFIG_ENABLE_BENCHMARKING'
  ]
)
test('tests',
     run_tests,
     workdir: meson.project_source_root(),
     # args: ['aod vs superliminal'],
     args: ['[aod::rtree]', '-d', 'yes']
     # args: ['[aod::rtree][fixme]']
     # args: ['drtree3 test: 5x5']
     # args: ['drtree3 test: 7x7']
     # args: ['drtree3 test: 10x10']
     # args: ['[basic test]']
     # args: ['[benchmark]']
     # args: ['[benchmark][aod::rtree]']
     # args: ['[drtree3]']
     # args: ['[drtree2]']
     # args: ['[benchmark]']
     # args: ['200x2d drtree']
     # args: ['allocations']
     # args: ['rtree basic ops']
     # args: ['bounds after removals & insertions']
    )

 

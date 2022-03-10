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
    RTree_dep,
  ],
  cpp_args: [
  ]
)
test('tests',
     run_tests,
     workdir: meson.project_source_root())


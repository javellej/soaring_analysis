project( 'soaring_analyser', 'cpp', default_options: ['cpp_std=c++1z', 'default_library=static']  )

compiler = meson.get_compiler('cpp')
cpp_lib = compiler.find_library( 'cpp-lib', dirs : meson.current_source_dir() + '/external/cpp-lib' )

#add_global_arguments( ['-O3', '-g', '-Wall', '-Werror', '-Wextra', '-fno-rtti'], language: 'cpp' )
add_global_arguments( ['-O3', '-g', '-Wall', '-Wextra', '-fno-rtti'], language: 'cpp' )
add_global_link_arguments( '-lpthread', language: 'cpp' )

includes = include_directories('external/cpp-lib/include','/usr/include/eigen3','/usr/include/png++')

executable( 'tracker', 'main.cpp',
            include_directories: [includes],
            dependencies: [cpp_lib],
            link_with: [] )

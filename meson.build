project('khpi_math',
    [ 'c', 'cpp' ],
    version: '1',
    default_options : ['cpp_std=c++20'])

cpp_compiler = meson.get_compiler('cpp')
cmake = import('cmake')

khpi_deps = []
khpi_incdirs = []

if host_machine.system() != 'emscripten'
    # vulkan
    khpi_deps += dependency('vulkan')
else
    # gl
    khpi_deps += dependency('gl')
endif

# glfw
if host_machine.system() != 'emscripten'
    khpi_deps += dependency('glfw3', static: true)
else
    # glfw_opts = cmake.subproject_options()
    # # glfw_opts.append_compile_args('cpp', '-DIMGUI_DISABLE_FILE_FUNCTIONS')
    # glfw_opts.append_link_args('-s USE_GLFW=3', '-s WASM=1', '-s ALLOW_MEMORY_GROWTH=1', '-s NO_EXIT_RUNTIME=0', '-s ASSERTIONS=1')
    # # glfw_opts.append_link_args('-lglfw3', language: 'cpp')
    # glfw_opts.add_cmake_defines({ 'EMSCRIPTEN': 'ON' })
    # glfw_proj = cmake.subproject('emscripten-glfw', options: glfw_opts)
    # glfw_dep = glfw_proj.dependency('glfw3')

    # glfw_override_dep = declare_dependency(
    #     include_directories: glfw_proj.include_directories('glfw3'),
    #     link_with: glfw_proj.target('glfw3')
    # )

    # meson.override_dependency('glfw3', glfw_override_dep)
    khpi_deps += dependency('sdl2', default_options: [ 'use_threads=disabled', 'use_atomic=disabled', 'assertions=disabled', 'test=false' ], static: true)
endif

# webgpu
# khpi_deps += dependency('wgpu-native')

# imgui
imgui_opts = [ 'default_library=static' ]
if host_machine.system() != 'emscripten'
    imgui_opts += [ 'glfw=enabled', 'vulkan=enabled' ]
else
    imgui_opts += [ 'sdl2=enabled', 'opengl=enabled' ]
endif

imgui_proj = subproject('imgui', default_options: imgui_opts)
khpi_deps += imgui_proj.get_variable('imgui_dep')

# implot
implot_proj = subproject('implot', default_options: [ 'default_library=static' ])
khpi_deps += implot_proj.get_variable('implot_dep')

khpi_sources = files(
    'src/main_window.cpp'
)

if host_machine.system() != 'emscripten'
    khpi_sources += files('src/main_vulkan.cpp')
else
    khpi_sources += files('src/main_webgpu.cpp')
endif

executable('khpi_math',
    khpi_sources,
    dependencies: khpi_deps,
    include_directories: khpi_incdirs
)
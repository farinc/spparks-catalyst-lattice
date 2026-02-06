py_mod = import('python')
py = py_mod.find_installation(required: true)

# Optional machine-style overrides
extra_cxxargs = get_option('extra_cxxargs')
extra_linkargs = get_option('extra_linkargs')
add_project_arguments(extra_cxxargs, language: 'cpp')
add_project_link_arguments(extra_linkargs, language: 'cpp')

# ---- Generate style_*.h (replacement for Make.sh style) ----
style_gen = custom_target(
  'gen_style_headers',
  input: files('../tools/meson/gen_styles.py'),
  output: [
    'style_app.h',
    'style_command.h',
    'style_diag.h',
    'style_dump.h',
    'style_pair.h',
    'style_region.h',
    'style_solve.h',
  ],
  command: [
    py, '@INPUT@',
    meson.current_source_dir(),
    meson.current_build_dir(),
  ],
  build_by_default: true,
)

inc_src = include_directories('.')

# ---- MPI / STUBS selection ----
mpi_feature = get_option('mpi')
mpi_dep = dependency(
  'mpi',
  language: 'cpp',
  required: (mpi_feature.enabled() or mpi_feature.auto())
)
use_mpi = mpi_dep.found()

stub_inc = include_directories('STUBS')
stub_lib = static_library(
  'spparks_mpi_stubs',
  files('STUBS/mpi.cpp'),
  include_directories: stub_inc,
  pic: true,
  build_by_default: false,
)

common_includes = inc_src
common_deps = []
common_link_with = []

if use_mpi
  common_deps += mpi_dep
else
  common_includes += stub_inc
  common_link_with += stub_lib
endif

# ---- OpenMP ----
openmp_feature = get_option('openmp')
openmp_dep = dependency('openmp', required: (openmp_feature.enabled() or openmp_feature.auto()))
if openmp_dep.found()
  common_deps += openmp_dep
endif

# ---- Hardcoded library sources (from src/Makefile.lib SRC=...) ----
lib_sources = files(
    'am_bezier.h',
    'am_ellipsoid.h',
    'am_raster.h',
    'am_teardrop.h',
    'app_am_ellipsoid.h',
    'app_chemistry.h',
    'app_diffusion.h',
    'app_diffusion_multiphase.h',
    'app_erbium.h',
    'app.h',
    'app_ising.h',
    'app_ising_single.h',
    'app_lattice.h',
    'app_membrane.h',
    'app_off_lattice.h',
    'app_phasefield_potts.h',
    'app_potts_am_bezier.h',
    'app_potts_am_path_gen.h',
    'app_potts_am_weld.h',
    'app_potts_grad.h',
    'app_potts.h',
    'app_potts_neigh.h',
    'app_potts_neighonly.h',
    'app_potts_pin.h',
    'app_potts_quaternion.h',
    'app_potts_strain.h',
    'app_potts_strain_pin.h',
    'app_potts_weld.h',
    'app_potts_weld_jom.h',
    'app_relax.h',
    'app_sinter.h',
    'app_sos.h',
    'app_test_group.h',
    'cluster.h',
    'comm_lattice.h',
    'comm_off_lattice.h',
    'create_box.h',
    'create_sites.h',
    'diag_array.h',
    'diag_cluster.h',
    'diag_diffusion.h',
    'diag_energy.h',
    'diag_erbium.h',
    'diag.h',
    'diag_propensity.h',
    'diag_sinter_density.h',
    'diag_sinter_free_energy.h',
    'diag_sinter_free_energy_pore.h',
    'domain.h',
    'dump.h',
    'dump_image.h',
    'dump_sites.h',
    'dump_text.h',
    'dump_vtk.h',
    'error.h',
    'finish.h',
    'fourth_order_bezier.h',
    'groups.h',
    'image.h',
    'input.h',
    'irregular.h',
    'kdtree.h',
    'lattice.h',
    'library.h',
    'math_const.h',
    'math_extra.h',
    'memory.h',
    'output.h',
    'pair.h',
    'pair_lj_cut.h',
    'pointers.h',
    'pool_shape.h',
    'potential.h',
    'potts_am_path_parser.h',
    'random_mars.h',
    'random_park.h',
    'read_sites.h',
    'region_block.h',
    'region_cylinder.h',
    'region.h',
    'region_intersect.h',
    'region_sphere.h',
    'region_union.h',
    'set.h',
    'shell.h',
    'solve_group.h',
    'solve.h',
    'solve_linear.h',
    'solve_tree.h',
    'spktype.h',
    'spparks.h',
    'teardrop.h',
    'timer.h',
    'universe.h',
    'variable.h',
    'variable_haz.h',
    'version.h',
    'weld_geometry.h',
    'app_surface_rxn.h',
    'am_bezier.cpp',
    'app_am_ellipsoid.cpp',
    'app_chemistry.cpp',
    'app.cpp',
    'app_diffusion.cpp',
    'app_diffusion_multiphase.cpp',
    'app_erbium.cpp',
    'app_ising.cpp',
    'app_ising_single.cpp',
    'app_lattice.cpp',
    'app_membrane.cpp',
    'app_off_lattice.cpp',
    'app_phasefield_potts.cpp',
    'app_potts_am_bezier.cpp',
    'app_potts_am_path_gen.cpp',
    'app_potts_am_weld.cpp',
    'app_potts.cpp',
    'app_potts_grad.cpp',
    'app_potts_neigh.cpp',
    'app_potts_neighonly.cpp',
    'app_potts_pin.cpp',
    'app_potts_quaternion.cpp',
    'app_potts_strain.cpp',
    'app_potts_strain_pin.cpp',
    'app_potts_weld.cpp',
    'app_potts_weld_jom.cpp',
    'app_relax.cpp',
    'app_sinter.cpp',
    'app_sos.cpp',
    'app_test_group.cpp',
    'cluster.cpp',
    'comm_lattice.cpp',
    'comm_off_lattice.cpp',
    'create_box.cpp',
    'create_sites.cpp',
    'diag_array.cpp',
    'diag_cluster.cpp',
    'diag.cpp',
    'diag_diffusion.cpp',
    'diag_energy.cpp',
    'diag_erbium.cpp',
    'diag_propensity.cpp',
    'diag_sinter_density.cpp',
    'diag_sinter_free_energy.cpp',
    'diag_sinter_free_energy_pore.cpp',
    'domain.cpp',
    'dump.cpp',
    'dump_image.cpp',
    'dump_sites.cpp',
    'dump_text.cpp',
    'dump_vtk.cpp',
    'error.cpp',
    'finish.cpp',
    'fourth_order_bezier.cpp',
    'groups.cpp',
    'image.cpp',
    'input.cpp',
    'irregular.cpp',
    'lattice.cpp',
    'library.cpp',
    'main.cpp',
    'math_extra.cpp',
    'memory.cpp',
    'output.cpp',
    'pair.cpp',
    'pair_lj_cut.cpp',
    'potential.cpp',
    'potts_am_path_parser.cpp',
    'random_mars.cpp',
    'random_park.cpp',
    'read_sites.cpp',
    'region_block.cpp',
    'region.cpp',
    'region_cylinder.cpp',
    'region_intersect.cpp',
    'region_sphere.cpp',
    'region_union.cpp',
    'set.cpp',
    'shell.cpp',
    'solve.cpp',
    'solve_group.cpp',
    'solve_linear.cpp',
    'solve_tree.cpp',
    'spparks.cpp',
    'timer.cpp',
    'universe.cpp',
    'variable.cpp',
    'app_surface_rxn.cpp'
) + [style_gen]

# ---- Libraries ----
libspparks_static = static_library(
  'spparks',
  lib_sources,
  include_directories: common_includes,
  dependencies: common_deps,
  link_with: common_link_with,
  install: get_option('build_static'),
)

libspparks_shared = shared_library(
  'spparks',
  lib_sources,
  include_directories: common_includes,
  dependencies: common_deps,
  link_with: common_link_with,
  install: get_option('build_shared'),
  soversion: '0',
)

# ---- Executable ----
if get_option('build_exe')
  executable(
    'spparks',
    files('main.cpp'),
    include_directories: common_includes,
    dependencies: common_deps,
    link_with: [libspparks_static] + common_link_with,
    install: true,
  )
endif

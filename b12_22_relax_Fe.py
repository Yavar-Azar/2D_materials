# -*- coding: utf-8 -*-
# -------------------------------------------------------------
# Bulk Configuration
# -------------------------------------------------------------

# Set up lattice
lattice = Triclinic(10.14*Angstrom, 5.84791*Angstrom, 7.09442*Angstrom, 90.0*Degrees, 90.0*Degrees, 90.0*Degrees)

# Define elements
elements = [Boron, Boron, Boron, Boron, Boron, Boron, Boron, Boron, Boron,
            Boron, Boron, Boron, Boron, Boron, Boron, Boron, Boron, Boron,
            Boron, Boron, Iron]

# Define coordinates
fractional_coordinates = [[ 0.51065163527 ,  0.003974884521, -0.000557267426],
                          [ 0.68005163527 ,  0.003974884521, -0.000557267426],
                          [ 0.84946163527 ,  0.003974884521, -0.000557267426],
                          [ 0.59536163527 ,  0.249414884521, -0.000557267426],
                          [ 0.76426163527 ,  0.249414884521, -0.000557267426],
                          [ 0.51065163527 ,  0.499924884521, -0.000557267426],
                          [ 0.68005163527 ,  0.499924884521, -0.000557267426],
                          [ 0.84946163527 ,  0.499924884521, -0.000557267426],
                          [ 0.59536163527 ,  0.745364884521, -0.000557267426],
                          [ 0.76426163527 ,  0.745364884521, -0.000557267426],
                          [ 0.00244163527 ,  0.003974884521, -0.000557267426],
                          [ 0.17184163527 ,  0.003974884521, -0.000557267426],
                          [ 0.34125163527 ,  0.003974884521, -0.000557267426],
                          [ 0.08714163527 ,  0.249414884521, -0.000557267426],
                          [ 0.25604163527 ,  0.249414884521, -0.000557267426],
                          [ 0.00244163527 ,  0.499924884521, -0.000557267426],
                          [ 0.17184163527 ,  0.499924884521, -0.000557267426],
                          [ 0.34125163527 ,  0.499924884521, -0.000557267426],
                          [ 0.08714163527 ,  0.745364884521, -0.000557267426],
                          [ 0.25604163527 ,  0.745364884521, -0.000557267426],
                          [ 0.420230529701,  0.266937840305,  0.144869675509]]

# Set up configuration
bulk_configuration = BulkConfiguration(
    bravais_lattice=lattice,
    elements=elements,
    fractional_coordinates=fractional_coordinates
    )

# -------------------------------------------------------------
# Calculator
# -------------------------------------------------------------
#----------------------------------------
# Basis Set
#----------------------------------------

# Basis set for Boron

boron_2p = ConfinedOrbital(
    principal_quantum_number=2,
    angular_momentum=1,
    radial_cutoff_radius=4.398046875*Angstrom,
    confinement_start_radius=3.5184375*Angstrom,
    additional_charge=0,
    confinement_strength=12.5*Hartree,
    confinement_power=2,
    radial_step_size=0.001*Bohr,
    )

boron_2p_0 = HydrogenOrbital(
    principal_quantum_number=2,
    angular_momentum=1,
    radial_cutoff_radius=4.9296875*Angstrom,
    confinement_start_radius=3.94375*Angstrom,
    charge=1.4,
    confinement_strength=12.5*Hartree,
    confinement_power=2,
    radial_step_size=0.001*Bohr,
    )

boron_2s = ConfinedOrbital(
    principal_quantum_number=2,
    angular_momentum=0,
    radial_cutoff_radius=4.52734375*Angstrom,
    confinement_start_radius=3.621875*Angstrom,
    additional_charge=0,
    confinement_strength=12.5*Hartree,
    confinement_power=2,
    radial_step_size=0.001*Bohr,
    )

boron_2s_0 = HydrogenOrbital(
    principal_quantum_number=2,
    angular_momentum=0,
    radial_cutoff_radius=3.580078125*Angstrom,
    confinement_start_radius=2.8640625*Angstrom,
    charge=4,
    confinement_strength=12.5*Hartree,
    confinement_power=2,
    radial_step_size=0.001*Bohr,
    )

boron_3d = HydrogenOrbital(
    principal_quantum_number=3,
    angular_momentum=2,
    radial_cutoff_radius=3.15234375*Angstrom,
    confinement_start_radius=2.521875*Angstrom,
    charge=4.8,
    confinement_strength=12.5*Hartree,
    confinement_power=2,
    radial_step_size=0.001*Bohr,
    )

BoronBasis = BasisSet(
    element=PeriodicTable.Boron,
    orbitals=[boron_2s, boron_2p, boron_2p_0, boron_3d, boron_2s_0],
    occupations=[5.0, -0.5, 0.5, -1.5, -0.5],
    hubbard_u=[0.0, 0.0, 0.0, 0.0, 0.0]*eV,
    dft_half_parameters=Automatic,
    filling_method=SphericalSymmetric,
    onsite_spin_orbit_split=[0.0, 0.0, 0.0, 0.0, 0.0]*eV,
    pseudopotential=NormConservingPseudoPotential("normconserving/pseudodojo/gga/standard/05_B.upf", local_potential_cutoff_threshold=1e-06*Hartree, local_potential_cutoff_radius=6.0*Bohr),
    )


# Basis set for Iron

iron_2p = HydrogenOrbital(
    principal_quantum_number=2,
    angular_momentum=1,
    radial_cutoff_radius=3.5068359375*Angstrom,
    confinement_start_radius=2.80546875*Angstrom,
    charge=2.2,
    confinement_strength=12.5*Hartree,
    confinement_power=2,
    radial_step_size=0.001*Bohr,
    )

iron_3d = ConfinedOrbital(
    principal_quantum_number=3,
    angular_momentum=2,
    radial_cutoff_radius=2.80859375*Angstrom,
    confinement_start_radius=2.246875*Angstrom,
    additional_charge=0,
    confinement_strength=12.5*Hartree,
    confinement_power=2,
    radial_step_size=0.001*Bohr,
    )

iron_3d_0 = HydrogenOrbital(
    principal_quantum_number=3,
    angular_momentum=2,
    radial_cutoff_radius=3.6091796875*Angstrom,
    confinement_start_radius=2.88734375*Angstrom,
    charge=3.1,
    confinement_strength=12.5*Hartree,
    confinement_power=2,
    radial_step_size=0.001*Bohr,
    )

iron_3p = ConfinedOrbital(
    principal_quantum_number=3,
    angular_momentum=1,
    radial_cutoff_radius=1.59765625*Angstrom,
    confinement_start_radius=1.278125*Angstrom,
    additional_charge=0,
    confinement_strength=12.5*Hartree,
    confinement_power=2,
    radial_step_size=0.001*Bohr,
    )

iron_3s = ConfinedOrbital(
    principal_quantum_number=3,
    angular_momentum=0,
    radial_cutoff_radius=1.41796875*Angstrom,
    confinement_start_radius=1.134375*Angstrom,
    additional_charge=0,
    confinement_strength=12.5*Hartree,
    confinement_power=2,
    radial_step_size=0.001*Bohr,
    )

iron_4f = HydrogenOrbital(
    principal_quantum_number=4,
    angular_momentum=3,
    radial_cutoff_radius=2.6875*Angstrom,
    confinement_start_radius=2.15*Angstrom,
    charge=9.4,
    confinement_strength=12.5*Hartree,
    confinement_power=2,
    radial_step_size=0.001*Bohr,
    )

iron_4s = ConfinedOrbital(
    principal_quantum_number=4,
    angular_momentum=0,
    radial_cutoff_radius=4.4296875*Angstrom,
    confinement_start_radius=3.54375*Angstrom,
    additional_charge=0,
    confinement_strength=12.5*Hartree,
    confinement_power=2,
    radial_step_size=0.001*Bohr,
    )

IronBasis = BasisSet(
    element=PeriodicTable.Iron,
    orbitals=[iron_3s, iron_3p, iron_4s, iron_3d, iron_4f, iron_2p, iron_3d_0],
    occupations=[1.0, 9.0, 0.0, 3.0, 1.5, 1.5, 0.0],
    hubbard_u=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]*eV,
    dft_half_parameters=Automatic,
    filling_method=SphericalSymmetric,
    onsite_spin_orbit_split=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]*eV,
    pseudopotential=NormConservingPseudoPotential("normconserving/pseudodojo/gga/stringent/26_Fe.upf", local_potential_cutoff_threshold=1e-06*Hartree, local_potential_cutoff_radius=6.0*Bohr),
    )

basis_set = [
    BoronBasis,
    IronBasis,
    ]

#----------------------------------------
# Exchange-Correlation
#----------------------------------------
exchange_correlation = SGGA.PBE

density_mesh_cutoff = OptimizedFFTGridSampling(
    energy_cutoff=200.0*Hartree,
    maximum_average_prime_factor=5.0,
    maximum_grid_enlargement=0.1,
    )

k_point_sampling = KpointDensity(
    density_a=1.0*Angstrom,
    density_b=1.0*Angstrom,
    density_c=0.0*Angstrom,
    symmetries=[
        ([[ 1., 0., 0.],
          [ 0., 1., 0.],
          [ 0., 0., 1.]], [ 0., 0., 0.]),
        ([[-1., 0., 0.],
          [ 0.,-1., 0.],
          [ 0., 0.,-1.]], [ 0., 0., 0.]),
        ],
    force_timereversal=True,
    shift_to_gamma=[True, True, True],
    )
exx_grid_cutoff = OptimizedFFTGridSampling(
    energy_cutoff=100.0*Hartree,
    maximum_average_prime_factor=5.0,
    maximum_grid_enlargement=0.1,
    )

exact_exchange_parameters = ExactExchangeParameters(
    aux_basis_tolerance=0.001,
    number_of_waves=1024,
    maximum_wave_vector=50.0,
    integral_tolerance=0.0001,
    relative_screening_tolerance=0.01,
    )
compensation_charge_mesh_cutoff = OptimizedFFTGridSampling(
    energy_cutoff=200.0*Hartree,
    maximum_average_prime_factor=5.0,
    maximum_grid_enlargement=0.1,
    )

numerical_accuracy_parameters = NumericalAccuracyParameters(
    density_mesh_cutoff=density_mesh_cutoff,
    k_point_sampling=k_point_sampling,
    radial_step_size=0.001*Bohr,
    density_cutoff=1e-06,
    interaction_max_range=20.0*Angstrom,
    number_of_reciprocal_points=1024,
    reciprocal_energy_cutoff=1250.0*Hartree,
    bands_per_electron=1.2,
    occupation_method=FermiDirac(1000.0*Kelvin*boltzmann_constant),
    exx_grid_cutoff=exx_grid_cutoff,
    exact_exchange_parameters=exact_exchange_parameters,
    compensation_charge_mesh_cutoff=compensation_charge_mesh_cutoff,
    )

iteration_control_parameters = IterationControlParameters(
    tolerance=1e-05,
    max_steps=100,
    algorithm=PulayMixer(),
    damping_factor=0.2,
    number_of_history_steps=5,
    start_mixing_after_step=0,
    mixing_variable=HamiltonianVariable,
    preconditioner=Kerker(energy_q0=0.01*Hartree, energy_qmax=200.0*Hartree, maximum_damping=0.01),
    linear_dependence_threshold=0.0,
    max_exx_updates=50,
    non_convergence_behavior=ContinueCalculation(),
    enable_scf_stop_file=True,
    )

poisson_solver = FastFourierSolver()

initialization_method = BasisSetInitialization()
density_matrix_method = PPCGSolver(
    absolute_tolerance=1e-08,
    relative_tolerance=1e-99,
    maximum_number_of_iterations=2,
    block_size=5,
    rr_period=3,
    buffer_states=0.05,
    initialization_method=initialization_method,
    )

algorithm_parameters = AlgorithmParameters(
    density_matrix_method=density_matrix_method,
    store_grids=True,
    store_basis_on_grid=Automatic,
    store_energy_density_matrix=Automatic,
    scf_restart_step_length=0.1*Angstrom,
    use_symmetries=Automatic,
    )

parallel_parameters = ParallelParameters(
    processes_per_neb_image=None,
    processes_per_individual=None,
    processes_per_bias_point=None,
    processes_per_saddle_search=1,
    )

checkpoint_handler = CheckpointHandler(
    time_interval=0.5*Hour,
    )

calculator = PlaneWaveCalculator(
    wave_function_cutoff=50.0*Hartree,
    basis_set=basis_set,
    exchange_correlation=exchange_correlation,
    numerical_accuracy_parameters=numerical_accuracy_parameters,
    iteration_control_parameters=iteration_control_parameters,
    poisson_solver=poisson_solver,
    algorithm_parameters=algorithm_parameters,
    parallel_parameters=parallel_parameters,
    checkpoint_handler=checkpoint_handler,
    charge=0.0,
    correction_extension=None,
    fixed_spin_moment=False,
    store_wave_functions=False,
    processes_per_kpoint=Automatic,
    )

bulk_configuration.setCalculator(calculator)
nlprint(bulk_configuration)
bulk_configuration.update()
nlsave('b12_22_vcrxed_Fe.hdf5', bulk_configuration)

# -------------------------------------------------------------
# Optimize Geometry
# -------------------------------------------------------------
bulk_configuration = OptimizeGeometry(
    bulk_configuration,
    max_forces=0.05*eV/Ang,
    max_steps=200,
    max_step_length=0.2*Ang,
    trajectory_filename='b12_22_vcrxed_Fe_trajectory.hdf5',
    trajectory_interval=5.0*Minute,
    restart_strategy=RestartFromTrajectory(),
    disable_stress=True,
    optimizer_method=LBFGS(),
    enable_optimization_stop_file=True,
)
nlsave('b12_22_vcrxed_Fe.hdf5', bulk_configuration)
nlprint(bulk_configuration)

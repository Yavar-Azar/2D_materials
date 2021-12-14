# -*- coding: utf-8 -*-
# -------------------------------------------------------------
# Bulk Configuration
# -------------------------------------------------------------

# Set up lattice
lattice = Triclinic(16.6858*Angstrom, 11.6807*Angstrom, 8.34289*Angstrom, 90.0*Degrees, 90.0*Degrees, 90.0*Degrees)

# Define elements
elements = [Boron, Boron, Boron, Boron, Boron, Boron, Boron, Boron, Boron,
            Boron, Boron, Boron, Boron, Boron, Boron, Boron, Boron, Boron,
            Boron, Boron, Boron, Boron, Boron, Boron, Boron, Boron, Boron,
            Boron, Boron, Boron, Boron, Boron, Boron, Boron, Boron, Boron,
            Boron, Boron, Boron, Boron, Boron, Boron, Boron, Boron, Boron,
            Boron, Boron, Boron, Boron, Boron, Boron, Boron, Boron, Boron,
            Boron, Boron, Boron, Boron, Boron, Boron, Boron, Boron, Boron,
            Boron, Iron]

# Define coordinates
fractional_coordinates = [[ 0.60463,  0.75277,  0.48781],
                          [ 0.70264,  0.75094,  0.50312],
                          [ 0.7988 ,  0.74907,  0.50002],
                          [ 0.89794,  0.74936,  0.50175],
                          [ 0.54804,  0.87435,  0.50207],
                          [ 0.649  ,  0.87655,  0.50134],
                          [ 0.85231,  0.8753 ,  0.50112],
                          [ 0.95189,  0.87546,  0.49898],
                          [ 0.60088,  0.0008 ,  0.50535],
                          [ 0.70159,  0.00063,  0.49247],
                          [ 0.79812,  0.     ,  0.49819],
                          [ 0.89743,  0.     ,  0.49878],
                          [ 0.5481 ,  0.12572,  0.50194],
                          [ 0.64891,  0.12394,  0.50233],
                          [ 0.85257,  0.1258 ,  0.50326],
                          [ 0.95186,  0.12632,  0.49754],
                          [ 0.60439,  0.24871,  0.49619],
                          [ 0.70271,  0.25079,  0.49785],
                          [ 0.79879,  0.24976,  0.50387],
                          [ 0.89785,  0.24977,  0.49881],
                          [ 0.54838,  0.37028,  0.50192],
                          [ 0.64993,  0.3733 ,  0.48663],
                          [ 0.85218,  0.37562,  0.49931],
                          [ 0.95223,  0.37611,  0.50115],
                          [ 0.60478,  0.50003,  0.51347],
                          [ 0.70284,  0.49991,  0.50139],
                          [ 0.79867,  0.49911,  0.49827],
                          [ 0.8985 ,  0.49981,  0.49901],
                          [ 0.54821,  0.63188,  0.49357],
                          [ 0.64952,  0.62662,  0.49503],
                          [ 0.85265,  0.62542,  0.49785],
                          [ 0.95212,  0.62498,  0.5033 ],
                          [ 0.10192,  0.74941,  0.50176],
                          [ 0.20113,  0.74904,  0.50004],
                          [ 0.29731,  0.75096,  0.5032 ],
                          [ 0.39535,  0.75275,  0.48749],
                          [ 0.04802,  0.87551,  0.49891],
                          [ 0.14758,  0.8753 ,  0.50112],
                          [ 0.35097,  0.87649,  0.50172],
                          [ 0.45193,  0.87436,  0.50187],
                          [ 0.10247, -0.     ,  0.49893],
                          [ 0.2018 , -0.     ,  0.49837],
                          [ 0.29833,  0.0006 ,  0.49284],
                          [ 0.39901,  0.00075,  0.50503],
                          [ 0.04804,  0.12628,  0.49756],
                          [ 0.14733,  0.12578,  0.50314],
                          [ 0.351  ,  0.12393,  0.50265],
                          [ 0.45181,  0.1257 ,  0.50172],
                          [ 0.10202,  0.24975,  0.49881],
                          [ 0.20109,  0.24985,  0.50376],
                          [ 0.29717,  0.25072,  0.49771],
                          [ 0.39552,  0.24868,  0.49591],
                          [ 0.04762,  0.37602,  0.50123],
                          [ 0.14774,  0.37561,  0.49933],
                          [ 0.34991,  0.37325,  0.4861 ],
                          [ 0.45157,  0.37027,  0.50107],
                          [ 0.1014 ,  0.4998 ,  0.49926],
                          [ 0.20119,  0.49903,  0.49833],
                          [ 0.29703,  0.4999 ,  0.50106],
                          [ 0.39517,  0.49995,  0.51093],
                          [ 0.04771,  0.62499,  0.5034 ],
                          [ 0.14717,  0.62545,  0.49776],
                          [ 0.35036,  0.62657,  0.49483],
                          [ 0.45169,  0.63197,  0.4927 ],
                          [ 0.49984,  0.50358,  0.6112 ]]

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
    density_a=2.0*Angstrom,
    density_b=2.0*Angstrom,
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
    occupation_method=MethfesselPaxton(1000.0*Kelvin*boltzmann_constant, 1),
    exx_grid_cutoff=exx_grid_cutoff,
    exact_exchange_parameters=exact_exchange_parameters,
    compensation_charge_mesh_cutoff=compensation_charge_mesh_cutoff,
    )

iteration_control_parameters = IterationControlParameters(
    tolerance=1e-05,
    max_steps=300,
    algorithm=PulayMixer(),
    damping_factor=0.2,
    number_of_history_steps=6,
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
    file_name='/home/yavar/MYGIT/2D_materials/borophene/x3_category/supercell/X3_23_Fesp.hdf5',
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
nlsave('X3_24Ferelaxed.hdf5', bulk_configuration)

# -------------------------------------------------------------
# Density Of States
# -------------------------------------------------------------
kpoint_grid = MonkhorstPackGrid(
    na=2,
    nb=3,
    )

density_of_states = DensityOfStates(
    configuration=bulk_configuration,
    kpoints=kpoint_grid,
    energy_zero_parameter=FermiLevel,
    bands_above_fermi_level=All,
    method=Full,
    )

nlsave('X3_24Ferelaxed.hdf5', density_of_states)
nlprint(density_of_states)

# -------------------------------------------------------------
# Electron Density
# -------------------------------------------------------------
electron_density = ElectronDensity(
    configuration=bulk_configuration,
    )
nlsave('/home/yavar/MYGIT/2D_materials/borophene/x3_category/supercell/chargefile.hdf5', electron_density)

# -------------------------------------------------------------
# Electrostatic Potential
# -------------------------------------------------------------
electrostatic_potential = ElectrostaticPotential(bulk_configuration)
nlsave('/home/yavar/MYGIT/2D_materials/borophene/x3_category/supercell/potentialfle.hdf5', electrostatic_potential)

# -------------------------------------------------------------
# Fat Bandstructure
# -------------------------------------------------------------
fat_bandstructure = FatBandstructure(
    configuration=bulk_configuration,
    route=['Y', 'G', 'X', 'V', 'U', 'R', 'Z', 'T'],
    points_per_segment=20,
    bands_above_fermi_level=All,
    projections=ProjectOnShellsByElement,
)
nlsave('/home/yavar/MYGIT/2D_materials/borophene/x3_category/supercell/fatbandfile.hdf5', fat_bandstructure)

# -------------------------------------------------------------
# Projected Density Of States
# -------------------------------------------------------------
kpoint_grid = MonkhorstPackGrid(
    na=2,
    nb=3,
    )

projected_density_of_states = ProjectedDensityOfStates(
    configuration=bulk_configuration,
    kpoints=kpoint_grid,
    projections=ProjectOnUpDownSpin,
    energies=numpy.linspace(-10, 10, 1001)*eV,
    energy_zero_parameter=FermiLevel,
    bands_above_fermi_level=All,
    spectrum_method=TetrahedronMethod,
)
nlsave('X3_24Ferelaxed.hdf5', projected_density_of_states)

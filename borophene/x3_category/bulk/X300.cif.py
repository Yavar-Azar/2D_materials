# Set up lattice
lattice = Triclinic(8.4*Angstrom, 2.9*Angstrom, 3.0*Angstrom, 90.0*Degrees, 90.0*Degrees, 90.0*Degrees)

# Define elements
elements = [Boron, Boron, Boron, Boron, Boron, Boron, Boron, Boron]

# Define coordinates
fractional_coordinates = [[ 0.2,  0. ,  0.5],
                          [ 0.4,  0. ,  0.5],
                          [ 0.6,  0. ,  0.5],
                          [ 0.8,  0. ,  0.5],
                          [ 0.1,  0.5,  0.5],
                          [ 0.3,  0.5,  0.5],
                          [ 0.7,  0.5,  0.5],
                          [ 0.9,  0.5,  0.5]]

# Set up configuration
bulk_configuration = BulkConfiguration(
    bravais_lattice=lattice,
    elements=elements,
    fractional_coordinates=fractional_coordinates
    )

OpenMM Gridded External Force Plugin
=====================

This project is a modified version of the [example plugin](https://github.com/openmm/openmmexampleplugin) plugin for [OpenMM](https://openmm.org).

This plugin defines a single Force subclass called GriddedExternalForce, which implements external forces based on cartesian grids. The energy and force equations therefore have the form $V(x,y,z)$ and $F(x,y,z)$, where the cartesian coordinates represent the position of target particles. Note that the plugin cannot construct these grids from forcefields, and so must be generated externally.

In practice, the GriddedExternalForce can be constructed in two ways:

1. With only the potential energy grid $V(x,y,z)$. Forces are then calculated using a second order center finite difference.
2. With the potential energy grid $V(x,y,z)$ and three force grids $F_x(x,y,z)$, $F_y(x,y,z)$, $F_z(x,y,z)$.

Values are always interpolated on the grid with a second order Taylor series expansion and first order finite difference derivative approximations (here the function $f$ either represents the potential $V$ or the force $F$):

$`
    f(\vec{r}) = f(\vec{r}_\text{grid}) + \Delta \vec{r}
    \begin{pmatrix}
    f^x\\
    f^y\\
    f^z\\
    \end{pmatrix}_{r_\text{grid}} + \frac{1}{2}\Delta \vec{r}
    \begin{pmatrix}
    f^{xx} & f^{xy} & f^{xz}\\
    f^{yx} & f^{yy} & f^{yz}\\
    f^{zx} & f^{zy} & f^{zz}\\
    \end{pmatrix}_{r_\text{grid}} \Delta \vec{r}^T 
`$

Building The Plugin
===================

Please refer to the [example plugin](https://github.com/openmm/openmmexampleplugin) instructions.

Python API
==========

The plugin is primarily built to be used with the Python API.

There are two ways to create this force:

Creating a GriddedExternalForce from NumPy arrays:

    from griddedexternalforceplugin import GriddedExternalForce
    
    # cap the force calculated by interpolation
    maxforce = 1e10 #kJ/mol/nm

    # create force
    force = GriddedExternalForce(Nx, Ny, Nz, V_grid, xmin, xmax, ymin, ymax, zmin, zmax, maxforce)

    # optionally, add the force grids
    force.setForcexGrid(Fx_grid)
    force.setForceyGrid(Fy_grid)
    force.setForcezGrid(Fz_grid)

    # add all particles to the force
    for i in range(nparticles):
        force.addParticle(i)

    # add force to system
    system = System()
    system.addForce(force)

This can be cumbersome to include in all scripts, and due to the way SWIG handles arrays, the script will require twice as much memory to store the grids (in Python and C++), which can be problematic for large grids.

Instead, the plugin offers the PotentialGridSetsFile Python object, which allows one to save their NumPy grids to a custom binary file format. The file path is then passed as an argument to create the GriddedExternalForce object, loading it directly in C++. 

It also allows for saving multiple sets of grids in one file, so that one file can contain the information of multiple GriddedExternalForce objects. For example, if we want to have an external force grid acting on the atoms of a water molecule, but want the H atoms to have different forces than the O atom, one would save these grids in the following way:

    from griddedexternalforceplugin import PotentialGridSetsFile

    # create the object
    vfile = PotentialGridSetsFile()

    # add the first set of grids
    vfile.add_grid_set(
        'Oxygen', # a string to label the grid set
        float_type=np.float64,
        Nx=xsize,
        Ny=ysize,
        Nz=zsize,
        xmin=xmin,
        xmax=xmax,
        ymin=ymin,
        ymax=ymax,
        zmin=zmin,
        zmax=zmax,
        V=V_grid_O,        
    )

    # add the second set of grids
    vfile.add_grid_set(
        'Hydrogen', # a string to label the grid set
        float_type=np.float64,
        Nx=xsize,
        Ny=ysize,
        Nz=zsize,
        xmin=xmin,
        xmax=xmax,
        ymin=ymin,
        ymax=ymax,
        zmin=zmin,
        zmax=zmax,
        V=V_grid_H,        
    )

    # write the grid sets to file
    vfile.write('potential_grid_sets_O_H.grd')

Now that we have saved this file, we can use it in any simulation to create GriddedExternalForce objects in the following way:

    from griddedexternalforceplugin import GriddedExternalForce, PotentialGridSetsFile

    verbose = True
    max_force = 1e10 #kJ/mol/nm

    gridfile_path = 'potential_grid_sets_O_H.grd'

    force_O = GriddedExternalForce(gridfile_path, 'Oxygen', max_force, verbose)
    force_H = GriddedExternalForce(gridfile_path, 'Hydrogen', max_force, verbose)

    # add the right particles to each force and add forces to the system...

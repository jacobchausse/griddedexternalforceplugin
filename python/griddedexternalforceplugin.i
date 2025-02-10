%module griddedexternalforceplugin

%import(module="openmm") "swig/OpenMMSwigHeaders.i"
%include "swig/typemaps.i"
%include "std_string.i"
/*
 * The following lines are needed to handle std::vector.
 * Similar lines may be needed for vectors of vectors or
 * for other STL types like maps.
 */

%include "std_vector.i"
namespace std {
  %template(vectord) vector<double>;
  %template(vectori) vector<int>;
};

%{
#include "GriddedExternalForce.h"
#include "OpenMM.h"
#include "OpenMMAmoeba.h"
#include "OpenMMDrude.h"
#include "openmm/RPMDIntegrator.h"
#include "openmm/RPMDMonteCarloBarostat.h"
%}

%pythoncode %{
import openmm as mm
import openmm.unit as unit
%}

/*
 * Convert C++ exceptions to Python exceptions.
*/
%exception {
    try {
        $action
    } catch (std::exception &e) {
        PyErr_SetString(PyExc_Exception, const_cast<char*>(e.what()));
        return NULL;
    }
}


namespace GriddedExternalForcePlugin {

class GriddedExternalForce : public OpenMM::Force {
public:
    GriddedExternalForce(int xsize, int ysize, int zsize, const std::vector<double>& potential, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double maxforce);

    GriddedExternalForce(std::string filename, std::string name, double maxforce, bool verbose);

    void GriddedExternalForce::setParameters(int xsize, int ysize, int zsize, const std::vector<double>& potential, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double maxforce);

    int GriddedExternalForce::addParticle(int particle);
    
    void GriddedExternalForce::setForcexGrid(const std::vector<double>& forcex);

    void GriddedExternalForce::setForceyGrid(const std::vector<double>& forcey);

    void GriddedExternalForce::setForcezGrid(const std::vector<double>& forcez);

    bool GriddedExternalForce::usesGriddedForce() const;

    /*
     * The reference parameters to this function are output values.
     * Marking them as such will cause swig to return a tuple.
    */
    %apply std::vector<int>& OUTPUT {std::vector<int>& particles};
    void GriddedExternalForce::getParticles(std::vector<int>& particles);
    %clear std::vector<int>& particles;    

    /*
     * Add methods for casting a Force to an GriddedExternalForce.
    */
    %extend {
        static GriddedExternalForcePlugin::GriddedExternalForce& cast(OpenMM::Force& force) {
            return dynamic_cast<GriddedExternalForcePlugin::GriddedExternalForce&>(force);
        }

        static bool isinstance(OpenMM::Force& force) {
            return (dynamic_cast<GriddedExternalForcePlugin::GriddedExternalForce*>(&force) != NULL);
        }
    }
};

}

%pythoncode %{
import numpy as np
from typing import Union, Optional, Literal, Dict

class PotentialGridSetsFile:
    
    _VERSION_NUMBER = 1

    _types_to_number = {np.float16: 0, np.float32: 1, np.float64: 2}
    _number_to_types = {0: np.float16, 1: np.float32, 2: np.float64}

    _types_to_str = {np.float16: 'e', np.float32: 'f', np.float64: 'd'}
    _str_to_type = {'e': np.float16, 'f': np.float32, 'd': np.float64}
    
    _endian_to_number = {'little': 0, 'big': 255}
    _number_to_endian = {0: 'little', 255: 'big'}
    
    _endian_to_str = {'little': '<', 'big': '>'}
    
    def __init__(self):
        """Create a new empty grid sets object.
        """

        self._ngridsets = 0
        self._names = []
        self._types = []
        self._V_list = []
        self._use_Fgrid_list = []
        self._Fx_list = []
        self._Fy_list = []
        self._Fz_list = []
        self._Nx_list = []
        self._Ny_list = []
        self._Nz_list = []
        self._xmin_list = []
        self._xmax_list = []
        self._ymin_list = []
        self._ymax_list = []
        self._zmin_list = []
        self._zmax_list = []

    def add_grid_set(self, 
        name: str,
        float_type: Union[np.float16, np.float32, np.float64],
        V: np.ndarray,
        Fx: Optional[np.ndarray],
        Fy: Optional[np.ndarray],
        Fz: Optional[np.ndarray],
        Nx: int,
        Ny: int,
        Nz: int,
        xmin: float,
        xmax: float,
        ymin: float,
        ymax: float,
        zmin: float,
        zmax: float
    ):
        """Add a grid set to the object.

        Parameters
        ----------
        name : str
            Name of the grid set. Needed to pick a grid set when loading in OpenMM.
        float_type : Union[np.float16, np.float32, np.float64]
            Float type to save the grid data in.
        V : np.ndarray
            Flattened grid for the potential energy (kJ/mol).
        Fx : Optional[np.ndarray]
            Flattened grid for the force in x (kJ/mol/nm).
        Fy : Optional[np.ndarray]
            Flattened grid for the force in y (kJ/mol/nm).
        Fz : Optional[np.ndarray]
            Flattened grid for the force in z (kJ/mol/nm).
        Nx : int
            Shape of the grid in x.
        Ny : int
            Shape of the grid in y.
        Nz : int
            Shape of the grid in z.
        xmin : float
            Lower limit position of the grid in x (nm).
        xmax : float
            Upper limit position of the grid in x (nm).
        ymin : float
            Lower limit position of the grid in y (nm).
        ymax : float
            Upper limit position of the grid in y (nm).
        zmin : float
            Lower limit position of the grid in z (nm).
        zmax : float
            Upper limit position of the grid in z (nm).
        """
        if float_type not in self._types_to_number.keys():
            raise Exception('Passed an incorrect type.')
        
        shape = V.shape
        if len(shape) != 1:
            raise Exception('V must be a 1-dimensional array.')

        use_F_grids = True
        if (Fx is None) and (Fx is None) and (Fx is None):
            use_F_grids = False

        if (Fx.shape != shape) and (Fy.shape != shape) and (Fz.shape != shape):
            raise Exception('Force grid shapes must match.')
        
        self._ngridsets += 1
        self._names.append(name)
        self._types.append(float_type)
        self._V_list.append(V)
        self._use_Fgrid_list.append(use_F_grids)
        self._Fx_list.append(Fx)
        self._Fy_list.append(Fy)
        self._Fz_list.append(Fz)
        self._Nx_list.append(Nx)
        self._Ny_list.append(Ny)
        self._Nz_list.append(Nz)
        self._xmin_list.append(xmin)
        self._xmax_list.append(xmax)
        self._ymin_list.append(ymin)
        self._ymax_list.append(ymax)
        self._zmin_list.append(zmin)
        self._zmax_list.append(zmax)

    def write(self, filepath: str, endian: Literal['big', 'little'] = 'little'):
        """Write the grid sets to a binary file (.grd)

        Parameters
        ----------
        filepath : str
            File path to save to.
        endian : Literal[\'big\', \'little\', optional
            Endianness of the saved data, by default 'little'. Only change if there are problems loading the data in OpenMM.
        """

        if endian not in ['big', 'little']:
            raise Exception('Endian must be either \'big\' or \'little\'.')

        if self._ngridsets == 0:
            raise Exception('No grid sets to write to file!')

        # Create the binary file
        with open(filepath, 'wb') as f:

            # Write endianness number
            f.write(self._endian_to_number[endian].to_bytes(1, endian))

            # Write version number
            f.write(self._VERSION_NUMBER.to_bytes(1, endian))

            # Number of grid sets
            f.write(self._ngridsets.to_bytes(1, endian))

            bytes_names = []
            length_header = 3 # endianness (uint8) + version number (uint8) + number of grid sets (uint8)

            for i, name in enumerate(self._names):
                b = bytearray()
                b.extend(name.encode('ASCII'))
                bytes_names.append(b)
                length_header += 1 # grid set number (uint8)
                length_header += 4 # length of grid set name (uint32)
                length_header += len(b) # grid set name (char8)
                length_header += 1 # (bool) if true, contains Force grids
                length_header += 1 # type of gris data (uint8). 0->float16, 1->float32, 2->float64
                length_header += 3*4 # grid shape (Nx, Ny, Nz) (uint32)
                length_header += 8 # position of grid set in file (uint64)
            
            position = int(length_header)

            for i in range(self._ngridsets):

                # Write grid set number
                f.write(i.to_bytes(1, endian))
                
                # Write length of name in bytes
                b = bytes_names[i]
                f.write(len(b).to_bytes(4, endian))

                # Write bytes of name
                f.write(b)

                # Write True is uses Force grids
                f.write(self._use_Fgrid_list[i].to_bytes(1, endian))

                # Write float type number
                f.write(self._types_to_number[self._types[i]].to_bytes(1, endian))

                # Write grid size
                Nx = self._Nx_list[i]
                Ny = self._Ny_list[i]
                Nz = self._Nz_list[i]
                f.write(Nx.to_bytes(4, endian))
                f.write(Ny.to_bytes(4, endian))
                f.write(Nz.to_bytes(4, endian))

                # Write position of the grid data in the file
                f.write(int(position).to_bytes(8, endian))

                type_size = np.dtype(self._types[i]).itemsize

                grids_byte_size = 0
                grids_byte_size += 6*type_size # x,y,z limits (depends on data type)
                grids_byte_size += Nx*Ny*Nz*type_size # potential grid data (depends on data type)

                if self._use_Fgrid_list[i]:
                    grids_byte_size += 3*Nx*Ny*Nz*type_size # x,y,z force grid data (depends on data type)

                position += grids_byte_size
            
            # Write grid sets
            for i in range(self._ngridsets):
                
                # Write x,y,z limits
                xmin = self._xmin_list[i]
                xmax = self._xmax_list[i]
                ymin = self._ymin_list[i]
                ymax = self._ymax_list[i]
                zmin = self._zmin_list[i]
                zmax = self._zmax_list[i]
                type_str = self._endian_to_str[endian] + self._types_to_str[self._types[i]]
                np.array([xmin, xmax, ymin, ymax, zmin, zmax]).astype(type_str).tofile(f)
                
                self._V_list[i].astype(type_str).tofile(f)

                if self._use_Fgrid_list[i]:
                    self._Fx_list[i].astype(type_str).tofile(f)
                    self._Fy_list[i].astype(type_str).tofile(f)
                    self._Fz_list[i].astype(type_str).tofile(f)

    @staticmethod
    def read_header(filepath: str) -> Dict:
        """Read only the header of a grid set binary file.

        Parameters
        ----------
        filepath : str
            Path to the file.

        Returns
        -------
        Dict
            Dictionary containing all the header information for the grid sets.
        """
        
        header_dict = {}
        
        # Open the binary file
        with open(filepath, 'rb') as f:

            # Read endianness number
            endian_number = int.from_bytes(f.read(1))
            if endian_number not in PotentialGridSetsFile._number_to_endian.keys():
                raise Exception(f'Could not read endianness (0 for little, 255 for big, got {endian_number}).')

            endian = PotentialGridSetsFile._number_to_endian[endian_number]

            header_dict['endian'] = endian

            # Read version number
            version_number = int.from_bytes(f.read(1), endian)
            header_dict['version'] = version_number

            # Read number of grid sets
            ngridsets = int.from_bytes(f.read(1), endian)
            header_dict['ngridsets'] = ngridsets
            
            header_dict['gridsets'] = []

            for i in range(ngridsets):
                grid_set_dict = {}

                # Read grid set number
                grid_set_number = int.from_bytes(f.read(1), endian)

                if i != grid_set_number:
                    raise Exception(f'Unmatched grid number in header (expected {i}, got {grid_set_number}).')
                
                grid_set_dict['number'] = grid_set_number

                # Read length of name in bytes
                name_length = int.from_bytes(f.read(4), endian)
                grid_set_dict['name_length'] = name_length

                # Read bytes of name
                name = f.read(name_length).decode('ASCII')
                grid_set_dict['name'] = name

                # Read uses Force grids
                uses_Fgrids = int.from_bytes(f.read(1), endian)
                if uses_Fgrids not in [0, 1]:
                    raise Exception(f'Got unexpected value for \'uses F grids\' ({uses_Fgrids}).')

                uses_Fgrids = bool(uses_Fgrids)
                grid_set_dict['uses_F_grids'] = uses_Fgrids

                # Read float type number
                type_number = int.from_bytes(f.read(1), endian)

                if type_number not in PotentialGridSetsFile._number_to_types.keys():
                    raise Exception(f'Unexpected value for type number ({type_number})')

                float_type = PotentialGridSetsFile._number_to_types[type_number]

                grid_set_dict['type'] = float_type

                # Read grid size
                Nx = int.from_bytes(f.read(4), endian)
                Ny = int.from_bytes(f.read(4), endian)
                Nz = int.from_bytes(f.read(4), endian)

                grid_set_dict['Nx'] = Nx
                grid_set_dict['Ny'] = Ny
                grid_set_dict['Nz'] = Nz

                position = int.from_bytes(f.read(8), endian)

                grid_set_dict['position'] = position

                header_dict['gridsets'].append(grid_set_dict)

        return header_dict
    
    @classmethod
    def read(cls, filepath: str):
        """Read a grid set file and store data in a new PotentialGridsFile object.

        Parameters
        ----------
        filepath : str
            Path to the grid set file.

        Returns
        -------
        PotentialGridsFile
            The created grid set object.
        """
        
        header_dict = PotentialGridSetsFile.read_header(filepath)
        pgridsfile = cls()
        endian = header_dict['endian']

        ngridsets = header_dict['ngridsets']

        # Open the binary file

        for i in range(ngridsets):

            with open(filepath, 'rb') as f:

                grid_set_dict = header_dict['gridsets'][i]

                float_type = grid_set_dict['type']
                shape = (grid_set_dict['Nx'], grid_set_dict['Ny'], grid_set_dict['Nz'])
                position = grid_set_dict['position']
                f.seek(position)
                
                type_size = np.dtype(grid_set_dict['type']).itemsize
                
                type_str = cls._endian_to_str[endian] + cls._types_to_str[grid_set_dict['type']]

                limits = np.fromfile(
                    f, 
                    dtype = type_str, 
                    count = 6
                )

                xmin = limits[0]
                xmax = limits[1]
                ymin = limits[2]
                ymax = limits[3]
                zmin = limits[4]
                zmax = limits[5]

                V = np.fromfile(
                    f, 
                    dtype = type_str, 
                    count = shape[0]*shape[1]*shape[2]
                )

                if grid_set_dict['uses_F_grids']:

                    Fx = np.fromfile(
                        f, 
                        dtype = type_str, 
                        count = shape[0]*shape[1]*shape[2]
                    )

                    Fy = np.fromfile(
                        f, 
                        dtype = type_str, 
                        count = shape[0]*shape[1]*shape[2]
                    )
                    
                    Fz = np.fromfile(
                        f, 
                        dtype = type_str, 
                        count = shape[0]*shape[1]*shape[2]
                    )
                
                else:
                    Fx = None
                    Fy = None
                    Fz = None

                pgridsfile.add_grid_set(
                    name = grid_set_dict['name'],
                    float_type = float_type,
                    V = V,
                    Fx = Fx,
                    Fy = Fy,
                    Fz = Fz,
                    Nx = shape[0],
                    Ny = shape[1],
                    Nz = shape[2],
                    xmin = xmin,
                    xmax = xmax,
                    ymin = ymin,
                    ymax = ymax,
                    zmin = zmin,
                    zmax = zmax
                )
        
        return pgridsfile
%}

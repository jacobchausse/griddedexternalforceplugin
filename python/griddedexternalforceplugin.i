%module griddedexternalforceplugin

%import(module="openmm") "swig/OpenMMSwigHeaders.i"
%include "swig/typemaps.i"

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

    void GriddedExternalForce::setParameters(int xsize, int ysize, int zsize, const std::vector<double>& potential, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double maxforce);

    int GriddedExternalForce::addParticle(int particle);
    
    void GriddedExternalForce::setForceGrids(const std::vector<double>& forcex, const std::vector<double>& forcey, const std::vector<double>& forcez);

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

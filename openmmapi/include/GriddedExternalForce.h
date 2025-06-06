#ifndef OPENMM_EXAMPLEFORCE_H_
#define OPENMM_EXAMPLEFORCE_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2014 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "openmm/Context.h"
#include "openmm/Force.h"
#include <vector>
#include "internal/windowsExportExample.h"

namespace GriddedExternalForcePlugin {

class OPENMM_EXPORT_EXAMPLE GriddedExternalForce : public OpenMM::Force {
public:
    /**
     * Create a gridded external force.
     * 
     * @param xsize     Shape of the grid in the x direction.
     * @param ysize     Shape of the grid in the y direction.
     * @param zsize     Shape of the grid in the z direction.
     * @param potential Flattened grid for the potential energy (kJ/mol).
     * @param xmin      Lower limit position of the grid in x (nm).
     * @param xmax      Upper limit position of the grid in x (nm).
     * @param ymin      Lower limit position of the grid in y (nm).
     * @param ymax      Upper limit position of the grid in y (nm).
     * @param zmin      Lower limit position of the grid in z (nm).
     * @param zmax      Upper limit position of the grid in z (nm).
     * @param maxforce  Force cutoff, useful to prevent force blow up in interpolation (kJ/mol/nm).
     */
    GriddedExternalForce(int xsize, int ysize, int zsize, const std::vector<double>& potential, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double maxforce);

    GriddedExternalForce(std::string filepath, std::string name, double maxforce, bool verbose=false);

    void setParameters(int xsize, int ysize, int zsize, const std::vector<double>& potential, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double maxforce);

    void getParameters(int& xsize, int& ysize, int& zsize, double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax, double& maxforce) const;

    void setForcexGrid(const std::vector<double>& forcex);

    void setForceyGrid(const std::vector<double>& forcey);

    void setForcezGrid(const std::vector<double>& forcez);

    void getGridPointers(const std::vector<double>*& ptrpotential, const std::vector<double>*& ptrforcex, const std::vector<double>*& ptrforcey, const std::vector<double>*& ptrforcez) const;

    int addParticle(int particle);

    bool usesGriddedForce() const;

    void getParticles(std::vector<int>& particles) const;

protected:
    OpenMM::ForceImpl* createImpl() const;

private:
    int xsize, ysize, zsize;
    double xmin, xmax, ymin, ymax, zmin, zmax, maxforce;
    std::vector<double> potential, forcex, forcey, forcez;
    std::vector<int> particles;
    bool griddedforcex, griddedforcey, griddedforcez;
};

} // namespace OpenMM


#endif /*OPENMM_EXAMPLEFORCE_H_*/

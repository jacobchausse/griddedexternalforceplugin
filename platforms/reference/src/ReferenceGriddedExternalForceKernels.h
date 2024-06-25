#ifndef REFERENCE_EXAMPLE_KERNELS_H_
#define REFERENCE_EXAMPLE_KERNELS_H_

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

#include "GriddedExternalForceKernels.h"
#include "openmm/Platform.h"
#include "openmm/Vec3.h"
#include <vector>

using namespace GriddedExternalForcePlugin;
using namespace OpenMM;
using namespace std;

/**
 * This kernel is invoked by GriddedExternalForce to calculate the forces acting on the system and the energy of the system.
 */
class ReferenceCalcGriddedExternalForceKernel : public CalcGriddedExternalForceKernel {
public:
    ReferenceCalcGriddedExternalForceKernel(string name, const Platform& platform) : CalcGriddedExternalForceKernel(name, platform) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the GriddedExternalForce this kernel will be used for
     */
    void initialize(const System& system, const GriddedExternalForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);

private:

    void taylorInterpolationDerivative(int i, int j, int k, double dxi, double dyi, double dzi, double& Fxi, double& Fyi, double& Fzi, double& Vi);

    void taylorInterpolation(const vector<double>& data, int i, int j, int k, double ic0, double jc0, double kc0, double& w);

    int xsize, ysize, zsize;
    double xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, dz, maxforce;
    const vector<double>* ptrpotential;
    const vector<double>* ptrforcex;
    const vector<double>* ptrforcey;
    const vector<double>* ptrforcez;
    vector<int> particles;
    bool griddedforce;
};

#endif /*REFERENCE_EXAMPLE_KERNELS_H_*/

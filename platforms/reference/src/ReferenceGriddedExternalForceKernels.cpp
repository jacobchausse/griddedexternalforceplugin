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

#include "ReferenceGriddedExternalForceKernels.h"
#include "GriddedExternalForce.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/Vec3.h"
#include "openmm/reference/ReferencePlatform.h"
#include <cmath>
#include <vector>

using namespace GriddedExternalForcePlugin;
using namespace OpenMM;
using namespace std;

static vector<Vec3>& extractPositions(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<Vec3>*) data->positions);
}

static vector<Vec3>& extractForces(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<Vec3>*) data->forces);
}

void ReferenceCalcGriddedExternalForceKernel::initialize(const System& system, const GriddedExternalForce& force) {    
    force.getParameters(xsize, ysize, zsize, xmin, xmax, ymin, ymax, zmin, zmax, maxforce);
    force.getParticles(particles);

    this->griddedforce = force.usesGriddedForce();
    force.getGridPointers(ptrpotential, ptrforcex, ptrforcey, ptrforcez);
    force.getPeriodicParameters(periodicx, periodicy, periodicz);

    xlen = (xmax - xmin);
    ylen = (ymax - ymin);
    zlen = (zmax - zmin);

    dx = xlen / (xsize-1);
    dy = ylen / (ysize-1);
    dz = zlen / (zsize-1);
}

double ReferenceCalcGriddedExternalForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<Vec3>& pos = extractPositions(context);
    vector<Vec3>& force = extractForces(context);
    int numParticles = particles.size();
    double energy = 0;
    
    // Compute the interactions.
    int i, j, k;
    double ic, jc, kc, ic0, jc0, kc0, Fx, Fy, Fz, V;
    double posx, posy, posz;
    for (int index = 0; index < numParticles; index++) {
        int particle = particles[index];
        
        posx = pos[particle][0];
        posy = pos[particle][1];
        posz = pos[particle][2];

        ic = (posx - xmin) / dx;
        jc = (posy - ymin) / dy;
        kc = (posz - zmin) / dz;

        i = ic + 0.5;
        j = jc + 0.5;
        k = kc + 0.5;
        
        // ensure index is within the boundary-1
        if (i < 1 || i > xsize-2 || j < 1 || j > ysize-2 || k < 1 || k > zsize-2)
            throw OpenMMException("GriddedExternalForce: position out of grid bounds.");

        ic0 = ic-i;
        jc0 = jc-j;
        kc0 = kc-k;

        if (griddedforce) {
            taylorInterpolation(*ptrpotential, i, j, k, ic0, jc0, kc0, V);
            taylorInterpolation(*ptrforcex, i, j, k, ic0, jc0, kc0, Fx);
            taylorInterpolation(*ptrforcey, i, j, k, ic0, jc0, kc0, Fy);
            taylorInterpolation(*ptrforcez, i, j, k, ic0, jc0, kc0, Fz);  
        }
        else {
            taylorInterpolationDerivative(i, j, k, ic0, jc0, kc0, Fx, Fy, Fz, V);
        }

        // cap the force
        Fx = max(-maxforce, min(Fx, maxforce));
        Fy = max(-maxforce, min(Fy, maxforce));
        Fz = max(-maxforce, min(Fz, maxforce));

        energy += V;

        force[particle] += Vec3(Fx, Fy, Fz);

    }
    return energy;
}

void ReferenceCalcGriddedExternalForceKernel::taylorInterpolation(const vector<double>& data, int i, int j, int k, double dxi, double dyj, double dzk, double& w) {
    
    // https://en.wikipedia.org/wiki/Finite_difference#Multivariate_finite_differences

    double dVdi0, dVdj0, dVdk0, d2Vdi20, d2Vdj20, d2Vdk20, d2Vdidj0, d2Vdidk0, d2Vdjdk0, Di, Dj, Dk;
    int ip1, im1, jp1, jm1, kp1, km1;

    if (periodicx) {
        ip1 = (i + 1) % (xsize - 1);
        im1 = (i - 1) % (xsize - 1);
        i = i % (xsize - 1);
    }
    else {  
        ip1 = (i + 1);
        im1 = (i - 1);
    }

    if (periodicy) {
        jp1 = (j + 1) % (ysize - 1);
        jm1 = (j - 1) % (ysize - 1);
        j = j % (ysize - 1);
    }
    else {  
        jp1 = (j + 1);
        jm1 = (j - 1);
    }

    if (periodicz) {
        kp1 = (k + 1) % (zsize - 1);
        km1 = (k - 1) % (zsize - 1);
        k = k % (zsize - 1);
    }
    else {  
        kp1 = (k + 1);
        km1 = (k - 1);
    }

    double wijk = data[k+zsize*(j+ysize*i)];
    double wip1 = data[k+zsize*(j+ysize*ip1)];
    double wim1 = data[k+zsize*(j+ysize*im1)];
    double wjp1 = data[k+zsize*(jp1+ysize*i)];
    double wjm1 = data[k+zsize*(jm1+ysize*i)];
    double wkp1 = data[kp1+zsize*(j+ysize*i)];
    double wkm1 = data[km1+zsize*(j+ysize*i)];

    dVdi0 = (wip1 - wim1) / (2);
    dVdj0 = (wjp1 - wjm1) / (2);
    dVdk0 = (wkp1 - wkm1) / (2);

    d2Vdi20 = (wip1 - 2*wijk + wim1);
    d2Vdj20 = (wjp1 - 2*wijk + wjm1);
    d2Vdk20 = (wkp1 - 2*wijk + wkm1);
    
    d2Vdidj0 = (data[k+zsize*(jp1+ysize*ip1)] - data[k+zsize*(jm1+ysize*ip1)] - data[k+zsize*(jp1+ysize*im1)] + data[k+zsize*(jm1+ysize*im1)])/(4);
    d2Vdidk0 = (data[kp1+zsize*(j+ysize*ip1)] - data[km1+zsize*(j+ysize*ip1)] - data[kp1+zsize*(j+ysize*im1)] + data[km1+zsize*(j+ysize*im1)])/(4);
    d2Vdjdk0 = (data[kp1+zsize*(jp1+ysize*i)] - data[kp1+zsize*(jm1+ysize*i)] - data[km1+zsize*(jp1+ysize*i)] + data[km1+zsize*(jm1+ysize*i)])/(4);

    Di = dxi*d2Vdi20  + dyj*d2Vdidj0 + dzk*d2Vdidk0;
    Dj = dxi*d2Vdidj0 + dyj*d2Vdj20  + dzk*d2Vdjdk0;
    Dk = dxi*d2Vdidk0 + dyj*d2Vdjdk0 + dzk*d2Vdk20 ;

    w = wijk + dxi*dVdi0 + dyj*dVdj0 + dzk*dVdk0 + 0.5*(dxi*Di + dyj*Dj + dzk*Dk);
}


void ReferenceCalcGriddedExternalForceKernel::taylorInterpolationDerivative(int i, int j, int k, double dxi, double dyj, double dzk, double& Fx, double& Fy, double& Fz, double& Vxyz) {
    // https://en.wikipedia.org/wiki/Finite_difference#Multivariate_finite_differences

    double dVdi0, dVdj0, dVdk0, d2Vdi20, d2Vdj20, d2Vdk20, d2Vdidj0, d2Vdidk0, d2Vdjdk0, Di, Dj, Dk;
    int ip1, im1, jp1, jm1, kp1, km1;

    if (periodicx) {
        ip1 = (i + 1) % (xsize - 1);
        im1 = (i - 1) % (xsize - 1);
        i = i % (xsize - 1);
    }
    else {  
        ip1 = (i + 1);
        im1 = (i - 1);
    }

    if (periodicy) {
        jp1 = (j + 1) % (ysize - 1);
        jm1 = (j - 1) % (ysize - 1);
        j = j % (ysize - 1);
    }
    else {  
        jp1 = (j + 1);
        jm1 = (j - 1);
    }

    if (periodicz) {
        kp1 = (k + 1) % (zsize - 1);
        km1 = (k - 1) % (zsize - 1);
        k = k % (zsize - 1);
    }
    else {  
        kp1 = (k + 1);
        km1 = (k - 1);
    }

    double Vijk = (*ptrpotential)[k+zsize*(j+ysize*i)];
    double Vip1 = (*ptrpotential)[k+zsize*(j+ysize*ip1)];
    double Vim1 = (*ptrpotential)[k+zsize*(j+ysize*im1)];
    double Vjp1 = (*ptrpotential)[k+zsize*(jp1+ysize*i)];
    double Vjm1 = (*ptrpotential)[k+zsize*(jm1+ysize*i)];
    double Vkp1 = (*ptrpotential)[kp1+zsize*(j+ysize*i)];
    double Vkm1 = (*ptrpotential)[km1+zsize*(j+ysize*i)];

    dVdi0 = (Vip1 - Vim1) / (2);
    dVdj0 = (Vjp1 - Vjm1) / (2);
    dVdk0 = (Vkp1 - Vkm1) / (2);

    d2Vdi20 = (Vip1 - 2*Vijk + Vim1);
    d2Vdj20 = (Vjp1 - 2*Vijk + Vjm1);
    d2Vdk20 = (Vkp1 - 2*Vijk + Vkm1);
    
    d2Vdidj0 = ((*ptrpotential)[k+zsize*(jp1+ysize*ip1)] - (*ptrpotential)[k+zsize*(jm1+ysize*ip1)] - (*ptrpotential)[k+zsize*(jp1+ysize*im1)] + (*ptrpotential)[k+zsize*(jm1+ysize*im1)])/(4);
    d2Vdidk0 = ((*ptrpotential)[kp1+zsize*(j+ysize*ip1)] - (*ptrpotential)[km1+zsize*(j+ysize*ip1)] - (*ptrpotential)[kp1+zsize*(j+ysize*im1)] + (*ptrpotential)[km1+zsize*(j+ysize*im1)])/(4);
    d2Vdjdk0 = ((*ptrpotential)[kp1+zsize*(jp1+ysize*i)] - (*ptrpotential)[kp1+zsize*(jm1+ysize*i)] - (*ptrpotential)[km1+zsize*(jp1+ysize*i)] + (*ptrpotential)[km1+zsize*(jm1+ysize*i)])/(4);

    Di = dxi*d2Vdi20  + dyj*d2Vdidj0 + dzk*d2Vdidk0;
    Dj = dxi*d2Vdidj0 + dyj*d2Vdj20  + dzk*d2Vdjdk0;
    Dk = dxi*d2Vdidk0 + dyj*d2Vdjdk0 + dzk*d2Vdk20 ;

    Fx = -(dVdi0 + Di)/dx;
    Fy = -(dVdj0 + Dj)/dy;
    Fz = -(dVdk0 + Dk)/dz;

    Vxyz = Vijk + dxi*dVdi0 + dyj*dVdj0 + dzk*dVdk0 + 0.5*(dxi*Di + dyj*Dj + dzk*Dk);
}
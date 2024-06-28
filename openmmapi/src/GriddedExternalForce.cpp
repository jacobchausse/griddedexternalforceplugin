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

#include "GriddedExternalForce.h"
#include "internal/GriddedExternalForceImpl.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/OpenMMException.h"
#include <vector>

using namespace GriddedExternalForcePlugin;
using namespace OpenMM;
using namespace std;

GriddedExternalForce::GriddedExternalForce(int xsize, int ysize, int zsize, const vector<double>& potential, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double maxforce) {
    setParameters(xsize, ysize, zsize, potential, xmin, xmax, ymin, ymax, zmin, zmax, maxforce);
}

void GriddedExternalForce::setParameters(int xsize, int ysize, int zsize, const vector<double>& potential, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double maxforce) {
    if (xsize < 2 || ysize < 2 || zsize < 2)
        throw OpenMMException("GriddedExternalForce: must have at least two points along each axis");
    if (potential.size() != xsize*ysize*zsize)
        throw OpenMMException("GriddedExternalForce: incorrect number of values in the potential");
    if (xmax <= xmin)
        throw OpenMMException("GriddedExternalForce: xmax <= xmin.");
    if (ymax <= ymin)
        throw OpenMMException("GriddedExternalForce: ymax <= ymin.");
    if (zmax <= zmin)
        throw OpenMMException("GriddedExternalForce: zmax <= zmin.");
    
    this->griddedforcex = false;
    this->griddedforcey = false;
    this->griddedforcez = false;
    this->potential = potential;
    this->xsize = xsize;
    this->ysize = ysize;
    this->zsize = zsize;
    this->xmin = xmin;
    this->xmax = xmax;
    this->ymin = ymin;
    this->ymax = ymax;
    this->zmin = zmin;
    this->zmax = zmax;
    this->maxforce = maxforce;
}

void GriddedExternalForce::getParameters(int& xsize, int& ysize, int& zsize, double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax, double& maxforce) const {
    xsize = this->xsize;
    ysize = this->ysize;
    zsize = this->zsize;
    xmin = this->xmin;
    xmax = this->xmax;
    ymin = this->ymin;
    ymax = this->ymax;
    zmin = this->zmin;
    zmax = this->zmax;
    maxforce = this->maxforce;
}

void GriddedExternalForce::setForcexGrid(const std::vector<double>& forcex) {
    if (forcex.size() != xsize*ysize*zsize)
        throw OpenMMException("GriddedExternalForce: incorrect number of values in forcex");
    
    this->griddedforcex = true;
    this->forcex = forcex;
}

void GriddedExternalForce::setForceyGrid(const std::vector<double>& forcey) {
    if (forcex.size() != xsize*ysize*zsize)
        throw OpenMMException("GriddedExternalForce: incorrect number of values in forcex");
    
    this->griddedforcey = true;
    this->forcey = forcey;
}

void GriddedExternalForce::setForcezGrid(const std::vector<double>& forcez) {
    if (forcex.size() != xsize*ysize*zsize)
        throw OpenMMException("GriddedExternalForce: incorrect number of values in forcex");
    
    this->griddedforcez = true;
    this->forcez = forcez;
}

void GriddedExternalForce::getGridPointers(const std::vector<double>*& ptrpotential, const std::vector<double>*& ptrforcex, const std::vector<double>*& ptrforcey, const std::vector<double>*& ptrforcez) const {
    ptrpotential = &(this->potential);
    ptrforcex = &(this->forcex);
    ptrforcey = &(this->forcey);
    ptrforcez = &(this->forcez);
}

int GriddedExternalForce::addParticle(int particle) {
    particles.push_back(particle);
    return particles.size()-1;
}

bool GriddedExternalForce::usesGriddedForce() const {
    bool missingforcegrid = (not ((not griddedforcex) && (not griddedforcey) && (not griddedforcez)))
                        and (not ((    griddedforcex) && (    griddedforcey) && (    griddedforcez)));

    if (missingforcegrid) throw OpenMMException("Not all force grids were set!");

    return griddedforcex && griddedforcey && griddedforcez;
}

void GriddedExternalForce::getParticles(vector<int>& particles) const {
    particles = this->particles;
}

ForceImpl* GriddedExternalForce::createImpl() const {
    return new GriddedExternalForceImpl(*this);
}
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
#include <iostream>
#include <fstream>

using namespace GriddedExternalForcePlugin;
using namespace OpenMM;
using namespace std;

const int GRDFILEVERSION = 1;

GriddedExternalForce::GriddedExternalForce(int xsize, int ysize, int zsize, const vector<double>& potential, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double maxforce) {
    setParameters(xsize, ysize, zsize, potential, xmin, xmax, ymin, ymax, zmin, zmax, maxforce);
}

GriddedExternalForce::GriddedExternalForce(std::string filepath, std::string name, double maxforce, bool verbose) {
    //open file
    std::ifstream file(filepath, std::ios::binary);
	
    //make sure file exists
    if (!file) {
		throw OpenMMException("Unable to open " + filepath + ".");
	}

    //read endianness
    uint8_t endianNumber;
    file.read(reinterpret_cast<char *>(&endianNumber), sizeof(endianNumber));

    int number = endianNumber;

    //read version number
    uint8_t versionNumber;
    file.read(reinterpret_cast<char *>(&versionNumber), sizeof(versionNumber));

    if (versionNumber != GRDFILEVERSION){
        throw OpenMMException("Incorrect grid file version.");
    }
    
    //read the number of grid sets
    uint8_t nGridSets_uint8;
    file.read(reinterpret_cast<char *>(&nGridSets_uint8), sizeof(nGridSets_uint8));  

    int nGridSets = nGridSets_uint8;

    bool matchedName = false;
    bool usesForceGrids;
    int typeNumber, Nx, Ny, Nz, position;
    uint8_t typeNumber_uint8, gridSetNumber_uint8;
    uint32_t Nx_uint32, Ny_uint32, Nz_uint32, nameLength_uint32;
    uint64_t position_uint64;
    

    for (int index = 0; index < nGridSets; index++) {
        
        //read the grid set number
        file.read(reinterpret_cast<char *>(&gridSetNumber_uint8), sizeof(gridSetNumber_uint8));
        int gridSetNumber = gridSetNumber_uint8;

        if (gridSetNumber != index) {
            throw OpenMMException("Unmatched grid set number.");
        }
        
        //read the name length
        file.read(reinterpret_cast<char *>(&nameLength_uint32), sizeof(nameLength_uint32));
        int nameLength = nameLength_uint32;

        //read the grid set name
        std::vector<char> gridSetName_char;
        gridSetName_char.resize(nameLength);

        file.read(&gridSetName_char[0], nameLength);

        std::string gridSetName(gridSetName_char.begin(), gridSetName_char.end());

        //if grid name doesnt match the one we want, skip to the next grid set header
        if (gridSetName != name) {
            file.seekg(2 + 12 + 8, std::ios::cur);
            continue;
        }

        matchedName = true;

        //read uses force grids
        file.read(reinterpret_cast<char *>(&usesForceGrids), 1);
        
        //read type number
        file.read(reinterpret_cast<char *>(&typeNumber_uint8), sizeof(typeNumber_uint8));
        typeNumber = typeNumber_uint8;

        //read Nx, Ny, Nz
        file.read(reinterpret_cast<char *>(&Nx_uint32), sizeof(Nx_uint32));
        file.read(reinterpret_cast<char *>(&Ny_uint32), sizeof(Ny_uint32));
        file.read(reinterpret_cast<char *>(&Nz_uint32), sizeof(Nz_uint32));
        Nx = Nx_uint32;
        Ny = Ny_uint32;
        Nz = Nz_uint32;

        //read grid data position
        file.read(reinterpret_cast<char *>(&position_uint64), sizeof(position_uint64));
        position = position_uint64;
    }

    if (!matchedName) {
        throw OpenMMException("Did not find \"" + name + "\" in the grid sets");
    }

    //seek to position of desired grid
    file.seekg(position, std::ios::beg);
    
    //only doubles are currently supported
    if (typeNumber == 2) {

        //read x,y,z limits
        double xmin1, xmax1, ymin1, ymax1, zmin1, zmax1;
        
        file.read(reinterpret_cast<char *>(&xmin1), sizeof(xmin1));
        file.read(reinterpret_cast<char *>(&xmax1), sizeof(xmax1));
        file.read(reinterpret_cast<char *>(&ymin1), sizeof(ymin1));
        file.read(reinterpret_cast<char *>(&ymax1), sizeof(ymax1));
        file.read(reinterpret_cast<char *>(&zmin1), sizeof(zmin1));
        file.read(reinterpret_cast<char *>(&zmax1), sizeof(zmax1));

        int gridSize = Nx*Ny*Nz;

        //read V data
        std::vector<double> V_in;
        V_in.resize(gridSize);

        char V_in_char[sizeof(double)];

        for (int j = 0; j < gridSize; j++) {
            file.read(&V_in_char[0], sizeof(double));
            std::copy(V_in_char, V_in_char + sizeof(double), reinterpret_cast<char*>(&V_in[j]));
            
        }     
        
        //set parameters in object
        setParameters(Nx, Ny, Nz, V_in, xmin1, xmax1, ymin1, ymax1, zmin1, zmax1, maxforce);  

        if (usesForceGrids) {
            
            char F_in_char[sizeof(double)];

            //read Fx data
            std::vector<double> Fx_in;
            Fx_in.resize(gridSize);

            for (int j = 0; j < gridSize; j++) {
                file.read(&F_in_char[0], sizeof(double));
                std::copy(F_in_char, F_in_char + sizeof(double), reinterpret_cast<char*>(&Fx_in[j]));
                
            }
            
            //read Fy data
            std::vector<double> Fy_in;
            Fy_in.resize(gridSize);

            for (int j = 0; j < gridSize; j++) {
                file.read(&F_in_char[0], sizeof(double));
                std::copy(F_in_char, F_in_char + sizeof(double), reinterpret_cast<char*>(&Fy_in[j]));
                
            }
            
            //read Fz data
            std::vector<double> Fz_in;
            Fz_in.resize(gridSize);

            for (int j = 0; j < gridSize; j++) {
                file.read(&F_in_char[0], sizeof(double));
                std::copy(F_in_char, F_in_char + sizeof(double), reinterpret_cast<char*>(&Fz_in[j]));
                
            } 

            //set force grids in object
            setForcexGrid(Fx_in);
            setForceyGrid(Fy_in);
            setForcezGrid(Fz_in);
        }
    }
    else {
        throw OpenMMException("Only doubles (float64) are currently supported.");
    }

    if (verbose) {
        double dx = (xmax - xmin) / (Nx-1);
        double dy = (ymax - ymin) / (Ny-1);
        double dz = (zmax - zmin) / (Nz-1);

        std::cout << std::endl << "Loaded Grid Set \"" << name << "\"" << std::endl
        << "Shape: (" << Nx << ", " << Ny << ", " << Nz << ")" << std::endl
        << "x Limits (nm): [" << xmin << ", " << xmax << "]" << ", dx=" << dx << std::endl
        << "y Limits (nm): [" << ymin << ", " << ymax << "]" << ", dy=" << dy <<std::endl 
        << "z Limits (nm): [" << zmin << ", " << zmax << "]" << ", dz=" << dz <<std::endl;
        
        if (usesForceGrids) std::cout << "Force grids: ON" << std::endl;
        else std::cout << "Force grids: OFF" << std::endl;
    }
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
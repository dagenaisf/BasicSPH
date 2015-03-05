/*******************************************************************************
 *
 * BasicSPH particle-based fluid solver
 * Copyright (C) 2015 François Dagenais
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Description: Basic SPH fluid solver.
 *
 ******************************************************************************/

#ifndef BASICSPH_H
#define BASICSPH_H

#include "SPHParticle.h"
#include "Vec3.h"

#include <vector>

class BasicSPH
{
public:
    BasicSPH(Vec3d volumeMin, Vec3d volumeMax, double mass, double restDensity, double h, double k, double dt);
    ~BasicSPH();

    std::vector<SPHParticle>& particles() { return _particles; }

    void enableAdaptiveTimeStep(double tolerance = 0.25, double min = 1E-7, double max = 1.0 );
    void disableAdaptativeTimeStep() { _useAdaptiveTimeStep = false; }
    void run(double time);

private:
    // Simulation steps
    void buildNeighbors();
    void computeDensityAndPressure();
    void addExternalForces();
    void computeArtificialViscosityForces();
    void computePressureForces();
    void integrate();

    // Adaptative time step
    void computeTimeStep();

    // Kernel (M4 - Cubic spline)
    void precomputeKernelCoefficients();
    double getKernelValue(double dist) const;
    Vec3d getKernelGradient(double dist, const Vec3d& xij) const;
    double getKernelLaplacian(double dist) const;

private:
    // Private structures
    struct Neighbor
    {
        Neighbor(long id, const Vec3d& xij, double dist) : id(id), dist(dist), xij(xij) {}
        long	id;
        double	dist;
        Vec3d	xij;
    };

    // Particles
    std::vector<SPHParticle>			_particles;
    std::vector<std::vector<Neighbor> >	_neighbors;

    // Simulation properties
    Vec3d	_volumeMin;		// Simulation volume min
    Vec3d	_volumeMax;		// Simulation volume max
    double	_mass;			// Mass of a single particle
    double	_restDensity;	// Rest density of the fluid
    double	_k;				// Pressure constant (gaz constant)
    double	_dt;			// Timestep size (in seconds)
    double	_h;				// Smoothing radius
    double	_bulkViscosity;	// Bulk viscosity
    double	_shearViscosity; // Shear viscosity

    // Adaptative time stepping
    bool	_useAdaptiveTimeStep;
    double	_tolerance;
    double	_minDT;
    double	_maxDT;
    double	_maxuij;

    // Kernel precomputed values & coefficients
    double	_halfH;
    double	_kernelValueCoeff;
    double	_kernelGradientCoeffA;
    double	_kernelGradientCoeffB;
    double	_kernelLaplacianCoeffA;
    double	_kernelLaplacianCoeffB;
};


#endif // BASICSPH_H

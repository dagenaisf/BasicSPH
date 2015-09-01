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
 * Description:	Basic SPH fluid solver based on the paper "Weakly compressible
 * 				SPH for free surface flows" (Becker and Teshner 2007). Surface
 * 				tension forces are not handled in this version, but they will
 * 				be in the future. We use a M4 Spline kernel in our calculations.
 * 				The particles are integrated using the semi-implicit Euler
 * 				method and we handle boundaries by repositionning particles
 * 				inside the simulation domain and by removing the velocity
 * 				component that moves the particle outside the domain.
 * 				This way of handling boundary conditions leads to particles
 * 				clustering near the boundaries of the simulation. A better
 * 				approach would be to use boundary particles with repulsive
 * 				forces. The adaptive time step calculation is based on the
 * 				paper "Smoothed Particle Hydrodynamics" (Monaghan 1992)
 *
 ******************************************************************************/

#include "BasicSPH.h"
#include "SpatialGrid.h"

#include <iostream>

namespace
{
    // This is much faster than calling pow(val, exponent)
    inline double pow2(double val) { return val*val; }
    inline double pow3(double val) { return val*val*val; }
    inline double pow7(double val) { return val*val*val*val*val*val*val; }
}

//--------------------------------------------------------------------------------------------------
// Constructor / Destructor
//--------------------------------------------------------------------------------------------------
BasicSPH::BasicSPH(Vec3d volumeMin, Vec3d volumeMax, double mass, double restDensity, double h,
                   double k, double dt)
    : _volumeMin(volumeMin),
      _volumeMax(volumeMax),
      _mass(mass),
      _restDensity(restDensity),
      _k(k),
      _dt(dt),
      _h(h),
      _bulkViscosity(0.5),
      _shearViscosity(0.0),
      _useAdaptiveTimeStep(false),
      _maxuij(0.0)
{
    precomputeKernelCoefficients();
}

BasicSPH::~BasicSPH()
{}

//--------------------------------------------------------------------------------------------------
// Public functions
//--------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------
// enableAdaptiveTimeStep
void BasicSPH::enableAdaptiveTimeStep(double tolerance, double min, double max)
{
    _useAdaptiveTimeStep = true;
    _tolerance = tolerance;
    _minDT = min;
    _maxDT = max;
}

//--------------------------------------------------------------------------------------------------
// run
void BasicSPH::run(double time)
{
    double oldTimeStep = _dt;

    // Run simulation!
    double timeLeft = time;
    while (timeLeft > 0.0)
    {
        // Run simulation steps
        buildNeighbors();
        computeDensityAndPressure();
        addExternalForces();
        computeArtificialViscosityForces();
        computePressureForces();

        // Compute time step
        if (_useAdaptiveTimeStep)
            computeTimeStep();

        // Limit timestep to the time left in the simulation
        if (timeLeft < _dt)
        {
            _dt = timeLeft;
        }

        // Update particles
        integrate();

        // Update time
        timeLeft -= _dt;

        std::cout << "Substep done! With timestep = " << _dt << std::endl;
    }

    // Restore old time step
    _dt = oldTimeStep;
}

//--------------------------------------------------------------------------------------------------
// Simulation stepsrhoij
//--------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------
// buildNeighbors
void BasicSPH::buildNeighbors()
{
    // Reserve space and initialize neighbors' data
    _neighbors.clear();
    _neighbors.resize(_particles.size());

    // Init spatial grid
    Vec3d borders(_h*2.0, _h*2.0, _h*2.0);
    Vec3d gridMin = _volumeMin;
    gridMin -= borders;
    Vec3d gridMax = _volumeMax;
    gridMax += borders;
    SpatialGrid<long> grid(_h, gridMin, gridMax);

    // Insert particles into grid
    for (long p=0; p<_particles.size(); ++p)
    {
        grid.insert(p, _particles[p].pos);
    }

    // Use grid to retrieve neighbor particles
    double h2 = _h*_h;
    std::vector<long*> nearbyParticles;
    for (long p=0; p<_particles.size(); ++p)
    {
        const SPHParticle &particle = _particles[p];

        // Get nearby particles
        grid.getElements(particle.pos, _h, nearbyParticles);

        // Find particles that are within smoothing radius
        _neighbors[p].reserve(50);
        for (long i=0; i<nearbyParticles.size(); ++i)
        {
            long nID = *nearbyParticles[i];
            const SPHParticle &neighborParticle = _particles[nID];

            // Skip current particle
            if (nID==p)
                continue;

            Vec3d xij = particle.pos - neighborParticle.pos;

            // Check if distance is lower than smoothing radius
            double dist2 = xij.dot(xij);
            if (dist2 < h2)
            {
                // Yup! Add the particle to the neighbors list along with
                // some precomputed informations
                _neighbors[p].push_back(Neighbor(nID, xij, sqrt(dist2)));
            }
        }
    }
}

//--------------------------------------------------------------------------------------------------
// computeDensityAndPressure
void BasicSPH::computeDensityAndPressure()
{
    // Precompute Taits coefficient
    const double B = (_k * _restDensity) / 7.0;

    // Iterate particles
    for (long p=0; p<_particles.size(); ++p)
    {
        SPHParticle &particle = _particles[p];

        // Reinitialize particle properties
        particle.density = 0.0;
        particle.oneOverDensity = 0.0;
        particle.accel = Vec3d(0,0,0);
        particle.pressure = 0.0;
        particle.oneOverDensity = 0.0;

        // Add current particle's contribution
        particle.density = _mass * getKernelValue(0.0);

        // Compute density from neighbors contributions
        // rho_i = SUM_j (m_j * Wij)
        for (long n=0; n<_neighbors[p].size(); ++n)
        {
            const Neighbor &neighbor = _neighbors[p][n];

            // Add contribution
            particle.density += _mass * getKernelValue(neighbor.dist);
        }

        // Precompute 1/rho and compute pressure using Tait's equation:
        // p_i = B * ((rho/rho_rest)^7-1)
        if (particle.density != 0.0)
        {
            particle.oneOverDensity = 1.0/particle.density;
            particle.pressure = B * (pow7(particle.density/_restDensity) - 1.0);
        }
    }
}

//--------------------------------------------------------------------------------------------------
// addExternalForces
void BasicSPH::addExternalForces()
{
    // Add gravity
    for (long p=0; p<_particles.size(); ++p)
    {
        SPHParticle &particle = _particles[p];
        particle.accel += Vec3d(0.0, -9.81, 0.0);
    }
}

//--------------------------------------------------------------------------------------------------
// computeArtificialViscosityForces
void BasicSPH::computeArtificialViscosityForces()
{
    const double alpha = _bulkViscosity;	// Bulk viscosity
    const double beta = _shearViscosity;	// Shear viscosity

    _maxuij = 0.0;

    // Precompute coefficients
    const double speedOfSound = sqrt(_k);	// c
    const double h2 = _h*_h;

    // Compute artificial viscosity forces
    for (long p=0; p<_particles.size(); ++p)
    {
        SPHParticle &particle = _particles[p];

        // No need to compute current particle's contribution, since its gradient is null!

        // Get neighbors contributions
        for (long n=0; n<_neighbors[p].size(); ++n)
        {
            const Neighbor &neighbor = _neighbors[p][n];
            const SPHParticle &neighborParticle = _particles[neighbor.id];

            // Compute contribution (based on the paper of Monaghan (1992))
            // fv_i/rho_i = SUM_j(m_j * IIij * gradient(Wij))
            //         | (alpha * c * uij + beta * uij^2) / avgRho,	when vij.xij < 0
            // IIij = -|
            //         | 0,											otherwise
            // uij = h * (vij.xij) / (|xij|^2 + 0.01*h^2)
            // vij = vi - vj
            // xij = xi - xj
            // avgRho = 0.5 * (rho_i + rho_j)
            Vec3d vij = particle.vel - neighborParticle.vel;
            double vijxij = vij.dot(neighbor.xij);
            double dij = neighbor.dist;
            double uij = _h*vijxij / (dij*dij + 0.01*h2);
            if (uij < 0)
            {
                // Compute contribution
                double avgDensity = 0.5 * (particle.density + neighborParticle.density);
                double IIij = (alpha*uij*speedOfSound + beta*uij*uij) / avgDensity;
                Vec3d contribution = getKernelGradient(neighbor.dist, neighbor.xij);
                contribution *= IIij;
                contribution *= _mass;

                particle.accel += contribution;
            }

            // Update maxuij for adaptive time step calculations
            if (uij > _maxuij)
            {
                _maxuij = uij;
            }
        }
    }
}

//--------------------------------------------------------------------------------------------------
// computePressureForces
void BasicSPH::computePressureForces()
{
    // Compute pressure forces
    for (long p=0; p<_particles.size(); ++p)
    {
        SPHParticle &particle = _particles[p];

        // No need to compute current particle's contribution, since its gradient is null!

        // Get neighbors contributions
        // fp_i/rho_i = -SUM_j (m_j * (pi/rho_i^2 + pj/rho_j^2) * gradient(Wij)
        for (long n=0; n<_neighbors[p].size(); ++n)
        {
            const Neighbor &neighbor = _neighbors[p][n];
            const SPHParticle &neighborParticle = _particles[neighbor.id];

            // Compute contribution
            Vec3d contribution = getKernelGradient(neighbor.dist, neighbor.xij);
            contribution *= (particle.pressure * particle.oneOverDensity * particle.oneOverDensity) +
                    (neighborParticle.pressure * neighborParticle.oneOverDensity * neighborParticle.oneOverDensity);
            contribution *= -1.0;
            contribution *= _mass;

            particle.accel += contribution;
        }
    }
}

//--------------------------------------------------------------------------------------------------
// integrate
void BasicSPH::integrate()
{
    // Update particles velocity and position
    for (long p=0; p<_particles.size(); ++p)
    {
        SPHParticle &particle = _particles[p];

        // Update velocity and position using the semi-implicit Euler method:
        // v(t+dt) = v(t) + a(t)*dt
        // x(t+dt) = x(t) + v(t+1)*dt
        particle.vel += particle.accel * _dt;
        particle.pos += particle.vel * _dt;

        // Apply boundary conditions
        if (particle.pos.x < _volumeMin.x)
        {
            particle.pos.x = _volumeMin.x;
            particle.vel.x = 0.0;
        }
        else if (particle.pos.x > _volumeMax.x)
        {
            particle.pos.x = _volumeMax.x;
            particle.vel.x = 0.0;
        }
        if (particle.pos.y < _volumeMin.y)
        {
            particle.pos.y = _volumeMin.y;
            particle.vel.y = 0.0;
        }
        else if (particle.pos.y > _volumeMax.y)
        {
            particle.pos.y = _volumeMax.y;
            particle.vel.y = 0.0;
        }
        if (particle.pos.z < _volumeMin.z)
        {
            particle.pos.z = _volumeMin.z;
            particle.vel.z = 0.0;
        }
        else if (particle.pos.z > _volumeMax.z)
        {
            particle.pos.z = _volumeMax.z;
            particle.vel.z = 0.0;
        }
    }
}

//--------------------------------------------------------------------------------------------------
// computeTimeStep
void BasicSPH::computeTimeStep()
{
    // Find maximum acceleration
    double maxAccel2 = 0.0;
    for (long p=0; p<_particles.size(); ++p)
    {
        const Vec3d &accel = _particles[p].accel;

        // Test squared acceleration length
        double accelLength2 = accel.x*accel.x + accel.y*accel.y + accel.z*accel.z;
        if (accelLength2 > maxAccel2)
        {
            maxAccel2 = accelLength2;
        }
    }

    // Compute force
    double maxForce = _mass * sqrt(maxAccel2);	// f = mass * a

    // Compute timestep (Based on paper "Smoothed Particles Hydrodynamics" (Monaghan 1992))
    // dt = tolerance * min(dta, dtcv)
    // dta = min_i(h/|f_i|)
    // dtcv = min_i( h/(c + 0.6*(alpha * c + beta * max_j(uij))) )
    double alpha = _bulkViscosity;
    double beta = _shearViscosity;
    double speedOfSound = sqrt(_k);	// c
    double tf = sqrt(_h / maxForce);
    double tcv = _h / (speedOfSound + 0.6*(alpha*speedOfSound + beta*_maxuij));
    _dt = (tf < tcv) ? tf : tcv;
    _dt *= _tolerance;

    // Clamp time step
    if (_dt < _minDT) _dt = _minDT;
    if (_dt > _maxDT) _dt = _maxDT;
}

//--------------------------------------------------------------------------------------------------
// Kernel functions
//--------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------
// precomputeKernelCoefficients
void BasicSPH::precomputeKernelCoefficients()
{
    const double PI = 3.14159265359;

    _halfH = _h/2.0;	// In Monaghan2005, h=half of smoothing radius

    // Precompute value coefficient (Identical for part A and B)
    _kernelValueCoeff = 1.0 / (4.0*PI*pow(_halfH,3));

    // Precompute gradient coefficients
    _kernelGradientCoeffA = 3.0 / (4.0*PI*pow(_halfH,4));
    _kernelGradientCoeffB = -3.0 / (4.0*PI*pow(_halfH,4));

    // Precompute laplacian coefficients
    _kernelLaplacianCoeffA = -9.0 / (PI*pow(_halfH,5));
    _kernelLaplacianCoeffB = 3.0 / (PI*pow(_halfH,5));

}

//--------------------------------------------------------------------------------------------------
// getKernelValue
inline
double BasicSPH::getKernelValue(double dist) const
{
    double q = dist/_halfH;
    if (q<1.0)
    {
        return _kernelValueCoeff * ( pow3(2.0-q)-4*pow3(1.0-q) );
    }
    else
    {
        return _kernelValueCoeff * pow3(2.0-q);
    }
}

//--------------------------------------------------------------------------------------------------
// getKernelGradient
inline
Vec3d BasicSPH::getKernelGradient(double dist, const Vec3d& xij) const
{
    double q = dist/_halfH;
    Vec3d gradient = xij;
    if (q<= 0.0)
    {
        gradient = Vec3d(0,0,0);
    }
    else if (q<1.0)
    {
        gradient *= _kernelGradientCoeffA * (4.0 * pow2(1.0-q) - pow2(2.0-q)) / dist;
    }
    else
    {
        gradient *= (_kernelGradientCoeffB * pow2(2.0 - q)) / dist;
    }

    return gradient;

}

//--------------------------------------------------------------------------------------------------
// getKernelLaplacian
inline
double BasicSPH::getKernelLaplacian(double dist) const
{
    double q = dist/_halfH;
    if (q<=1.0)
    {
        return _kernelLaplacianCoeffA * (1.0-q);
    }
    else
    {
        return _kernelLaplacianCoeffB * (3.0-q-(2.0/q));
    }
}

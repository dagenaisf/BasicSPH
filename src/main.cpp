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
 * Description: Main class that uses the BasicSPH solver to simulate a fluid.
 * 				The result is output in .geo format after every frame.
 *
 ******************************************************************************/

#include "BasicSPH.h"
#include "SPHParticleHoudiniIO.h"

#include <string>
#include <math.h>
#include <iostream>

//--------------------------------------------------------------------------------------------------
// initParticlesBox
void initParticlesBox(const Vec3d& 				boxMin,
                      const Vec3d& 				boxMax,
                      double 					spacing,
                      std::vector<SPHParticle>&	particles)
{
    int resX = static_cast<int>(floor((boxMax.x-boxMin.x) / spacing));
    int resY = static_cast<int>(floor((boxMax.y-boxMin.y) / spacing));
    int resZ = static_cast<int>(floor((boxMax.z-boxMin.z) / spacing));

    // Clear previous particles
    particles.clear();

    // Fill box with particles
    particles.reserve(resX*resY*resZ);
    for (int x=0; x<resX; ++x)
    {
        for (int y=0; y<resY; ++y)
        {
            for (int z=0; z<resZ; ++z)
            {
                SPHParticle particle;

                // Init position
                particle.pos.x = static_cast<double>(x) * spacing + boxMin.x;
                particle.pos.y = static_cast<double>(y) * spacing + boxMin.y;
                particle.pos.z = static_cast<double>(z) * spacing + boxMin.z;

                // Init other properties
                particle.vel = Vec3d(0,0,0);

                // Add particle
                particles.push_back(particle);
            }
        }
    }
}

//--------------------------------------------------------------------------------------------------
// MAIN
int main(int argc, char* argv[])
{
    // Validate arguments
    if (argc != 2)
    {
        std::cout << "Usage: " << argv[0] << " OutputDirectory" << std::endl;
        return -1;
    }

    std::string outputPath(argv[1]);


    // Init simulator
    std::cout << "Initializing simulator..." << std::endl;
    const Vec3d volumeMin(-3,   0, -1);
    const Vec3d volumeMax( 3,   3,  1);
    const double mass = 1.5;
    const double restDensity = 998.23;
    const double h = 0.2;
    const double k = 100.0;
    const double dt = 0.001;
    BasicSPH sph(volumeMin, volumeMax, mass, restDensity, h, k, dt);
    sph.enableAdaptiveTimeStep();

    // Init particles
    std::cout << "Initializing particles" << std::endl;
    Vec3d boxMin(-3, 0, -1);
    Vec3d boxMax(-1, 2,  1);
    initParticlesBox(boxMin, boxMax, 0.1, sph.particles());

    // Output first frame (initial frames)
    const std::string filePrefix("particles_");
    SPHParticleHoudiniIO::outputParticles(sph.particles(), outputPath, filePrefix, 1);

    // Run simulation and output a frame every 1/24 second
    std::cout << "Running simulation!" << std::endl;
    const double frameTime = 1.0/24.0;
    const double totalSimulationTime = 10.0;
    double time = 0.0;
    int currentFrame = 2;
    while (time < totalSimulationTime)
    {
        std::cout << "Simulating frame " << currentFrame << " (" << (frameTime+time) << "s)";
        std::cout << std::endl;

        // Run simulation
        sph.run(frameTime);

        // Output particles
        SPHParticleHoudiniIO::outputParticles(sph.particles(), outputPath, filePrefix, currentFrame);

        // Update simulation time
        time += frameTime;
        ++currentFrame;
    }

    std::cout << "Done!" << std::endl;

    return 0;
}


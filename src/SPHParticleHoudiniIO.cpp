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
 * Description: Adapter class between BasicSPH and HoudiniFileDumpHelper
 *
 ******************************************************************************/

#include "SPHParticleHoudiniIO.h"
#include "HoudiniFileDumpHelper.h"
#include "SPHParticle.h"

#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>

namespace
{
    //----------------------------------------------------------------------------------------------
    // Data provider: This class is used to provide point data (attributes) to the
    //                HoudiniFilDumpHelper class.
    class SPHParticleDataProvider : public HoudiniFileDumpParticles::ParticlesDataProvider
    {
    public:
        SPHParticleDataProvider(const std::vector<SPHParticle>& particles)
            : _particles(particles)
        {}

        virtual int getNbPoints() { return _particles.size(); }
        virtual int getNbAttributes() { return 4; }

        virtual void getAttributesInfo(std::vector<std::string>	&names,
                                       std::vector<int>			&nbValues)
        {
            names.push_back("v");
            names.push_back("accel");
            names.push_back("density");
            names.push_back("pressure");

            nbValues.push_back(3);	// v
            nbValues.push_back(3);	// accel
            nbValues.push_back(1);	// density
            nbValues.push_back(1);	// pressure
        }

        virtual void getPtPosition(int ptID, float& posX, float& posY, float& posZ, float& posW)
        {
            const SPHParticle &particle = _particles[ptID];
            posX = particle.pos.x;
            posY = particle.pos.y;
            posZ = particle.pos.z;
            posW = 1.0;
        }

        virtual void getPtAttributes(int ptID, const std::string&	attribDelim,
                                     const std::string&				valueDelim,
                                     std::ostream&					out)
        {
            const SPHParticle& p = _particles[ptID];
            out << p.vel.x << valueDelim << p.vel.y << valueDelim << p.vel.z << attribDelim;
            out << p.accel.x << valueDelim << p.accel.y << valueDelim << p.accel.z << attribDelim;
            out << p.density << attribDelim;
            out << p.pressure << attribDelim;
        }

    private:
        const std::vector<SPHParticle> &_particles;
    };
}

//--------------------------------------------------------------------------------------------------
// Public functions
//--------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------
// outputParticles
void SPHParticleHoudiniIO::outputParticles(const std::vector<SPHParticle>&	particles,
                                           const std::string&				outputPath,
                                           const std::string&				prefix,
                                           int								frameNo)
{
    // Build output file pathname
    std::stringstream pathname;
    pathname << outputPath;
    if (*outputPath.rbegin() != '/' && *outputPath.rbegin() != '\\')
        pathname << "/";
    pathname << prefix << frameNo << ".geo";

    // Open file
    std::ofstream file;
    file.open(pathname.str().c_str(), std::ofstream::out | std::ofstream::trunc);
    if (!file.is_open())
    {
        std::cout << "Error: Failed to open file '" << pathname.str() << "' for writting..." << std::endl;
        return;
    }

    // Output particles
    SPHParticleDataProvider dataProvider(particles);
    HoudiniFileDumpParticles particlesDumper(&dataProvider);
    particlesDumper.dump(file);

    // Close file
    file.close();
}


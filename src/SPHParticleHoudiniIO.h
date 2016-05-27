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

#ifndef SPHPARTICLEHOUDINIIO_H
#define SPHPARTICLEHOUDINIIO_H

#include <vector>
#include <string>

struct SPHParticle;

namespace SPHParticleHoudiniIO
{
    void outputParticles(const std::vector<SPHParticle>&	particles,
                         const std::string&					outputPath,
                         const std::string&					prefix,
                         int								frameNo);
}

#endif // SPHPARTICLEHOUDINIIO_H

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
 * Description: SPH particle data structure
 *
 ******************************************************************************/

#ifndef SPHPARTICLE_H
#define SPHPARTICLE_H

#include "Vec3.h"

struct SPHParticle
{
    Vec3d	pos;
    Vec3d	vel;
    Vec3d	accel;
    double	density;
    double	oneOverDensity;
    double	pressure;
};

#endif // SPHPARTICLE_H

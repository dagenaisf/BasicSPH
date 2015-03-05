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
 * Description: Helper class used to output particles' data in the .geo format
 * 				that can be read by the Houdini software (http://www.sidefx.com).
 * 				This class is a generic implementation for any type of particles.
 *
 *
 ******************************************************************************/

#ifndef HOUDINIFILEDUMPHELPER_H
#define HOUDINIFILEDUMPHELPER_H

#include <string>
#include <vector>
#include <ostream>

class HoudiniFileDumpParticles
{
public:
	class ParticlesDataProvider;

public:
	HoudiniFileDumpParticles();
	HoudiniFileDumpParticles(ParticlesDataProvider* dataProvider);
	~HoudiniFileDumpParticles();

	void setDataProvider(ParticlesDataProvider* dataProvider);
	void dump(std::ostream& stream);

private:
	ParticlesDataProvider*	_dataProvider;
};

class HoudiniFileDumpParticles::ParticlesDataProvider
{
public:

	ParticlesDataProvider() {}
	virtual ~ParticlesDataProvider() {}

	virtual int getNbPoints() = 0;
	virtual int getNbAttributes() = 0;
	virtual void getAttributesInfo(std::vector<std::string>& names, std::vector<int>& nbValues) = 0;

	virtual void getPtPosition(int ptID, float& posX, float& posY, float& posZ, float& posW) = 0;
	virtual void getPtAttributes(int ptID, const std::string& attribDelim, const std::string& valueDelim, std::ostream& out) = 0;
};

#endif	// HOUDINIFILEDUMPHELPER_H

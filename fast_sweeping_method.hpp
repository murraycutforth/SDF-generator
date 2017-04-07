/*
 * DESCRIPTION
 * 	These functions implement the fast sweeping method in 3D to close the zero level set
 * 	surface, and set the values away from the narrow band to the appropriate sign.
 * 	They each consist of 8 iterations over every single cell, in various orderings such
 * 	that all characteristic directions are covered.
 * 
 * AUTHOR
 * 	Murray Cutforth
 * 
 * Date
 * 	07/04/2017
 */

#ifndef FSM_H

#include "APSS.hpp"
#include "Nvector.hpp"
#include <vector>

typedef std::vector<std::vector<std::vector<float> > > floatarray3D;
typedef std::vector<std::vector<std::vector<int> > > intarray3D;


void fast_sweeping_method_fillbulkcells (

	floatarray3D& sdfarray, 
	const intarray3D& activecells, 
	unsigned int Nx, 
	unsigned int Ny, 
	unsigned int Nz, 
	float dx, 
	float dy, 
	float dz, 
	float UNSETVAL, 
	int UBOUND
);

void fast_sweeping_method_closesurface (

	floatarray3D& sdfarray, 
	intarray3D& activecells, 
	const APSS& kenneth, 
	int Nx, 
	int Ny, 
	int Nz, 
	float dx, 
	float dy, 
	float dz, 
	float x0, 
	float y0, 
	float z0, 
	int UBOUND,
	float UNSETVAL
);

#endif

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

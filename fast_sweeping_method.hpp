#ifndef FSM_H

#include <vector>

typedef std::vector<std::vector<std::vector<float> > > floatarray3D;
typedef std::vector<std::vector<std::vector<int> > > intarray3D;


void fast_sweeping_method (floatarray3D& sdfarray, const intarray3D& activecells, unsigned int Nx, unsigned int Ny, unsigned int Nz, float dx, float dy, float dz, float UNSETVAL, int UBOUND);

#endif

#include "fast_sweeping_method.hpp"
#include <algorithm>
#include <cmath>
#include <cassert>
#include <vector>




bool comp (float a, float b)
{
	return fabs(a) < fabs(b);
}




void FSM_cellupdate (floatarray3D& sdfarray, int i, int j, int k, int Nx, int Ny, int Nz, float UNSETVAL)
{
	float a_m = sdfarray[std::max(i-1,0)][j][k];
	float a = sdfarray[std::min(i+1,Nx-1)][j][k];
	float b_m = sdfarray[i][std::max(0,j-1)][k];
	float b = sdfarray[i][std::min(Ny-1,j+1)][k];
	float c_m = sdfarray[i][j][std::max(k-1,0)];
	float c = sdfarray[i][j][std::min(k+1,Nz-1)];


	if (fabs(a_m) < fabs(a))
	{
		a = a_m;
		a_m = fabs(a);
	}
	else
	{
		a_m = fabs(a);
	}

	
	if (fabs(b_m) < fabs(b))
	{
		b = b_m;
		b_m = fabs(b);
	}
	else
	{
		b_m = fabs(b);
	}


	if (fabs(c_m) < fabs(c))
	{
		c = c_m;
		c_m = fabs(c);
	}
	else
	{
		c_m = fabs(c);
	}

	std::vector<float> sdfs (3, 0.0);
	sdfs.push_back(a);
	sdfs.push_back(b);
	sdfs.push_back(c);

	std::sort(sdfs.begin(), sdfs.end(), comp);
	sdfarray[i][j][k] = sdfs.front();
}


void fast_sweeping_method (floatarray3D& sdfarray, const intarray3D& activecells, unsigned int Nx, unsigned int Ny, unsigned int Nz, float dx, float dy, float dz, float UNSETVAL, int UBOUND)
{
	for (int i=0; i<int(Nx); ++i)
	{
		for (int j=0; j<int(Ny); ++j)
		{
			for (int k=0; k<int(Nz); ++k)
			{
				if (activecells[i][j][k] != UBOUND)
				{
					FSM_cellupdate(sdfarray,i,j,k,Nx,Ny,Nz,UNSETVAL);
				}
			}
		}
	}

	for (int i=0; i<int(Nx); ++i)
	{
		for (int j=int(Ny)-1; j>=0; --j)
		{
			for (int k=0; k<int(Nz); ++k)
			{
				if (activecells[i][j][k] != UBOUND)
				{
					FSM_cellupdate(sdfarray,i,j,k,Nx,Ny,Nz,UNSETVAL);
				}
			}
		}
	}

	for (int i=0; i<int(Nx); ++i)
	{
		for (int j=int(Ny)-1; j>=0; --j)
		{
			for (int k=int(Nz)-1; k>=0; --k)
			{
				if (activecells[i][j][k] != UBOUND)
				{
					FSM_cellupdate(sdfarray,i,j,k,Nx,Ny,Nz,UNSETVAL);
				}
			}
		}
	}

	for (int i=0; i<int(Nx); ++i)
	{
		for (int j=0; j<int(Ny); ++j)
		{
			for (int k=int(Nz)-1; k>=0; --k)
			{
				if (activecells[i][j][k] != UBOUND)
				{
					FSM_cellupdate(sdfarray,i,j,k,Nx,Ny,Nz,UNSETVAL);
				}
			}
		}
	}
	
	for (int i=int(Nx)-1; i>=0; --i)
	{
		for (int j=0; j<int(Ny); ++j)
		{
			for (int k=0; k<int(Nz); ++k)
			{
				if (activecells[i][j][k] != UBOUND)
				{
					FSM_cellupdate(sdfarray,i,j,k,Nx,Ny,Nz,UNSETVAL);
				}
			}
		}
	}

	for (int i=int(Nx)-1; i>=0; --i)
	{
		for (int j=int(Ny)-1; j>=0; --j)
		{
			for (int k=0; k<int(Nz); ++k)
			{
				if (activecells[i][j][k] != UBOUND)
				{
					FSM_cellupdate(sdfarray,i,j,k,Nx,Ny,Nz,UNSETVAL);
				}
			}
		}
	}

	for (int i=int(Nx)-1; i>=0; --i)
	{
		for (int j=int(Ny)-1; j>=0; --j)
		{
			for (int k=int(Nz)-1; k>=0; --k)
			{
				if (activecells[i][j][k] != UBOUND)
				{
					FSM_cellupdate(sdfarray,i,j,k,Nx,Ny,Nz,UNSETVAL);
				}
			}
		}
	}

	for (int i=int(Nx)-1; i>=0; --i)
	{
		for (int j=0; j<int(Ny); ++j)
		{
			for (int k=int(Nz)-1; k>=0; --k)
			{
				if (activecells[i][j][k] != UBOUND)
				{
					FSM_cellupdate(sdfarray,i,j,k,Nx,Ny,Nz,UNSETVAL);
				}
			}
		}
	}
}

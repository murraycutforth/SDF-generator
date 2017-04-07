#include "fast_sweeping_method.hpp"
#include <algorithm>
#include <cmath>
#include <cassert>
#include <iostream>
#include <vector>




bool comp (float a, float b)
{
	return fabs(a) < fabs(b);
}




void FSM_cellupdate (floatarray3D& sdfarray, const intarray3D& activecells, int i, int j, int k, int Nx, int Ny, int Nz, float UNSETVAL)
{
	
	int iP = std::min(Nx-1,i+1);
	int iM = std::max(0,i-1);
	int jP = std::min(Ny-1,j+1);
	int jM = std::max(0,j-1);
	int kP = std::min(Nz-1,k+1);
	int kM = std::max(0,k-1);
	
	std::vector<float> sdfs;
	
	if (activecells[iP][j][k] != 0 && sdfarray[iP][j][k] != UNSETVAL) sdfs.push_back(sdfarray[iP][j][k]);
	if (activecells[iM][j][k] != 0 && sdfarray[iM][j][k] != UNSETVAL) sdfs.push_back(sdfarray[iM][j][k]);
	if (activecells[i][jP][k] != 0 && sdfarray[i][jP][k] != UNSETVAL) sdfs.push_back(sdfarray[i][jP][k]);
	if (activecells[i][jM][k] != 0 && sdfarray[i][jM][k] != UNSETVAL) sdfs.push_back(sdfarray[i][jM][k]);
	if (activecells[i][j][kP] != 0 && sdfarray[i][j][kP] != UNSETVAL) sdfs.push_back(sdfarray[i][j][kP]);
	if (activecells[i][j][kM] != 0 && sdfarray[i][j][kM] != UNSETVAL) sdfs.push_back(sdfarray[i][j][kM]);
	
	if (sdfs.empty()) sdfs.push_back(sdfarray[iP][j][k]);
	
	float avgval = 0.0;
	for (float el : sdfs) avgval += el;
	avgval /= sdfs.size();
	sdfarray[i][j][k] = avgval;
}




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
)
{
	for (int i=0; i<int(Nx); ++i)
	{
		for (int j=0; j<int(Ny); ++j)
		{
			for (int k=0; k<int(Nz); ++k)
			{
				if (activecells[i][j][k] == UBOUND)
				{
					FSM_cellupdate(sdfarray,activecells,i,j,k,Nx,Ny,Nz,UNSETVAL);
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
				if (activecells[i][j][k] == UBOUND)
				{
					FSM_cellupdate(sdfarray,activecells,i,j,k,Nx,Ny,Nz,UNSETVAL);
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
				if (activecells[i][j][k] == UBOUND)
				{
					FSM_cellupdate(sdfarray,activecells,i,j,k,Nx,Ny,Nz,UNSETVAL);
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
				if (activecells[i][j][k] == UBOUND)
				{
					FSM_cellupdate(sdfarray,activecells,i,j,k,Nx,Ny,Nz,UNSETVAL);
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
				if (activecells[i][j][k] == UBOUND)
				{
					FSM_cellupdate(sdfarray,activecells,i,j,k,Nx,Ny,Nz,UNSETVAL);
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
				if (activecells[i][j][k] == UBOUND)
				{
					FSM_cellupdate(sdfarray,activecells,i,j,k,Nx,Ny,Nz,UNSETVAL);
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
				if (activecells[i][j][k] == UBOUND)
				{
					FSM_cellupdate(sdfarray,activecells,i,j,k,Nx,Ny,Nz,UNSETVAL);
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
				if (activecells[i][j][k] == UBOUND)
				{
					FSM_cellupdate(sdfarray,activecells,i,j,k,Nx,Ny,Nz,UNSETVAL);
				}
			}
		}
	}
}











void closesurface_cellupdate(floatarray3D& sdfarray, intarray3D& activecells, const APSS& kenneth, Nvector& pos, int i, int j, int k, int Nx, int Ny, int Nz, float dx, float dy, float dz, int UBOUND, float UNSETVAL)
{
	// Compare sign of sdfarray[i][j][k] with sign of surrounding 6 adjacent cells which also have status==1
	
	bool open = false;
	
	int iP = std::min(Nx-1,i+1);
	int iM = std::max(0,i-1);
	int jP = std::min(Ny-1,j+1);
	int jM = std::max(0,j-1);
	int kP = std::min(Nz-1,k+1);
	int kM = std::max(0,k-1);
	
	if (activecells[iP][j][k] == 1 && std::copysign(1.0,sdfarray[iP][j][k]) != std::copysign(1.0,sdfarray[i][j][k])) open = true;
	if (activecells[iM][j][k] == 1 && std::copysign(1.0,sdfarray[iM][j][k]) != std::copysign(1.0,sdfarray[i][j][k])) open = true;
	if (activecells[i][jP][k] == 1 && std::copysign(1.0,sdfarray[i][jP][k]) != std::copysign(1.0,sdfarray[i][j][k])) open = true;
	if (activecells[i][jM][k] == 1 && std::copysign(1.0,sdfarray[i][jM][k]) != std::copysign(1.0,sdfarray[i][j][k])) open = true;
	if (activecells[i][j][kP] == 1 && std::copysign(1.0,sdfarray[i][j][kP]) != std::copysign(1.0,sdfarray[i][j][k])) open = true;
	if (activecells[i][j][kM] == 1 && std::copysign(1.0,sdfarray[i][j][kM]) != std::copysign(1.0,sdfarray[i][j][k])) open = true;
	
	
	// If any of these signs are different, compute SDF in all surrounding cells with status == UBOUND
	
	if (open)
	{		
		float h = 2.0*fabs(sdfarray[i][j][k]);
		
		if (activecells[iP][j][k] == UBOUND)
		{
			pos.data[0] += dx;
			sdfarray[iP][j][k] = kenneth.evaluate_surface(pos,h);
			pos.data[0] -= dx;
			activecells[iP][j][k] = 1;
		}
		
		if (activecells[iM][j][k] == UBOUND)
		{
			pos.data[0] -= dx;
			sdfarray[iM][j][k] = kenneth.evaluate_surface(pos,h);
			pos.data[0] += dx;
			activecells[iM][j][k] = 1;
			
		}
		
		if (activecells[i][jP][k] == UBOUND)
		{
			pos.data[1] += dy;
			sdfarray[i][jP][k] = kenneth.evaluate_surface(pos,h);
			pos.data[1] -= dy;
			activecells[i][jP][k] = 1;
			
		}
		
		if (activecells[i][jM][k] == UBOUND) 
		{
			pos.data[1] -= dy;
			sdfarray[i][jM][k] = kenneth.evaluate_surface(pos,h);
			pos.data[1] += dy;
			activecells[i][jM][k] = 1;
			
		}
		
		if (activecells[i][j][kP] == UBOUND) 
		{
			pos.data[2] += dz;
			sdfarray[i][j][kP] = kenneth.evaluate_surface(pos,h);
			pos.data[2] -= dz;
			activecells[i][j][kP] = 1;
			
		}
		
		if (activecells[i][j][kM] == UBOUND)
		{
			pos.data[2] -= dz;
			sdfarray[i][j][kM] = kenneth.evaluate_surface(pos,h);
			pos.data[2] += dz;
			activecells[i][j][kM] = 1;
			
		}
	}
	
	activecells[i][j][k] = 0;
}




void fast_sweeping_method_closesurface (floatarray3D& sdfarray, intarray3D& activecells, const APSS& kenneth, int Nx, int Ny, int Nz, float dx, float dy, float dz, float x0, float y0, float z0, int UBOUND, float UNSETVAL)
{
	Nvector pos (3,0.0);
	float x,y,z;
	
	for (int i=0; i<Nx; ++i)
	{
		x = x0 + 0.5*dx + i*dx;
		
		for (int j=0; j<Ny; ++j)
		{
			y = y0 + 0.5*dy + j*dy;
			
			for (int k=0; k<Nz; ++k)
			{
				if (activecells[i][j][k] == 1)
				{
					if (sdfarray[i][j][k] == UNSETVAL) std::cout << "ERROR: unsetval in cell " << i << "," << j << "," << k << std::endl;
					assert(sdfarray[i][j][k] != UNSETVAL);
					
					z = z0 + 0.5*dz + k*dz;
					pos.data[0] = x;
					pos.data[1] = y;
					pos.data[2] = z;
					
					closesurface_cellupdate(sdfarray,activecells,kenneth,pos,i,j,k,Nx,Ny,Nz,dx,dy,dz,UBOUND,UNSETVAL);
				}
			}
		}
	}
	
	std::cout << "\rSweep 1 complete.." << std::flush;

	for (int i=0; i<Nx; ++i)
	{
		x = x0 + 0.5*dx + i*dx;
		
		for (int j=Ny-1; j>=0; --j)
		{
			y = y0 + 0.5*dy + j*dy;
			
			for (int k=0; k<Nz; ++k)
			{
				if (activecells[i][j][k] == 1)
				{
					assert(sdfarray[i][j][k] != UNSETVAL);
					
					z = z0 + 0.5*dz + k*dz;
					pos.data[0] = x;
					pos.data[1] = y;
					pos.data[2] = z;
					
					closesurface_cellupdate(sdfarray,activecells,kenneth,pos,i,j,k,Nx,Ny,Nz,dx,dy,dz,UBOUND,UNSETVAL);
				}
			}
		}
	}
	
	std::cout << "\rSweep 2 complete.." << std::flush;

	for (int i=0; i<Nx; ++i)
	{
		x = x0 + 0.5*dx + i*dx;
		
		for (int j=Ny-1; j>=0; --j)
		{
			y = y0 + 0.5*dy + j*dy;
			
			for (int k=Nz-1; k>=0; --k)
			{
				if (activecells[i][j][k] == 1)
				{
					assert(sdfarray[i][j][k] != UNSETVAL);
					
					z = z0 + 0.5*dz + k*dz;
					pos.data[0] = x;
					pos.data[1] = y;
					pos.data[2] = z;
					
					closesurface_cellupdate(sdfarray,activecells,kenneth,pos,i,j,k,Nx,Ny,Nz,dx,dy,dz,UBOUND,UNSETVAL);
				}
			}
		}
	}
	
	std::cout << "\rSweep 3 complete.." << std::flush;

	for (int i=0; i<Nx; ++i)
	{
		x = x0 + 0.5*dx + i*dx;
		
		for (int j=0; j<Ny; ++j)
		{
			y = y0 + 0.5*dy + j*dy;
			
			for (int k=Nz-1; k>=0; --k)
			{
				if (activecells[i][j][k] == 1)
				{
					assert(sdfarray[i][j][k] != UNSETVAL);
					
					z = z0 + 0.5*dz + k*dz;
					pos.data[0] = x;
					pos.data[1] = y;
					pos.data[2] = z;
					
					closesurface_cellupdate(sdfarray,activecells,kenneth,pos,i,j,k,Nx,Ny,Nz,dx,dy,dz,UBOUND,UNSETVAL);
				}
			}
		}
	}
	
	std::cout << "\rSweep 4 complete.." << std::flush;
	
	for (int i=Nx-1; i>=0; --i)
	{
		x = x0 + 0.5*dx + i*dx;
		
		for (int j=0; j<Ny; ++j)
		{
			y = y0 + 0.5*dy + j*dy;
			
			for (int k=0; k<Nz; ++k)
			{
				if (activecells[i][j][k] == 1)
				{
					assert(sdfarray[i][j][k] != UNSETVAL);
					
					z = z0 + 0.5*dz + k*dz;
					pos.data[0] = x;
					pos.data[1] = y;
					pos.data[2] = z;
					
					closesurface_cellupdate(sdfarray,activecells,kenneth,pos,i,j,k,Nx,Ny,Nz,dx,dy,dz,UBOUND,UNSETVAL);
				}
			}
		}
	}
	
	std::cout << "\rSweep 5 complete.." << std::flush;

	for (int i=Nx-1; i>=0; --i)
	{
		x = x0 + 0.5*dx + i*dx;
		
		for (int j=Ny-1; j>=0; --j)
		{
			y = y0 + 0.5*dy + j*dy;
			
			for (int k=0; k<Nz; ++k)
			{
				if (activecells[i][j][k] == 1)
				{
					assert(sdfarray[i][j][k] != UNSETVAL);
					
					z = z0 + 0.5*dz + k*dz;
					pos.data[0] = x;
					pos.data[1] = y;
					pos.data[2] = z;
					
					closesurface_cellupdate(sdfarray,activecells,kenneth,pos,i,j,k,Nx,Ny,Nz,dx,dy,dz,UBOUND,UNSETVAL);
				}
			}
		}
	}
	
	std::cout << "\rSweep 6 complete.." << std::flush;

	for (int i=Nx-1; i>=0; --i)
	{
		x = x0 + 0.5*dx + i*dx;
		
		for (int j=Ny-1; j>=0; --j)
		{
			y = y0 + 0.5*dy + j*dy;
			
			for (int k=Nz-1; k>=0; --k)
			{
				if (activecells[i][j][k] == 1)
				{
					assert(sdfarray[i][j][k] != UNSETVAL);
					
					z = z0 + 0.5*dz + k*dz;
					pos.data[0] = x;
					pos.data[1] = y;
					pos.data[2] = z;
					
					closesurface_cellupdate(sdfarray,activecells,kenneth,pos,i,j,k,Nx,Ny,Nz,dx,dy,dz,UBOUND,UNSETVAL);
				}
			}
		}
	}
	
	std::cout << "\rSweep 7 complete.." << std::flush;

	for (int i=Nx-1; i>=0; --i)
	{
		x = x0 + 0.5*dx + i*dx;
		
		for (int j=0; j<Ny; ++j)
		{
			y = y0 + 0.5*dy + j*dy;
			
			for (int k=Nz-1; k>=0; --k)
			{
				if (activecells[i][j][k] == 1)
				{
					assert(sdfarray[i][j][k] != UNSETVAL);
					
					z = z0 + 0.5*dz + k*dz;
					pos.data[0] = x;
					pos.data[1] = y;
					pos.data[2] = z;
					
					closesurface_cellupdate(sdfarray,activecells,kenneth,pos,i,j,k,Nx,Ny,Nz,dx,dy,dz,UBOUND,UNSETVAL);
				}
			}
		}
	}
	
	std::cout << "\rSweep 8 complete.." << std::flush;
	
	for (int i=Nx-1; i>=0; --i)
	{		
		for (int j=0; j<Ny; ++j)
		{
			for (int k=Nz-1; k>=0; --k)
			{
				if (activecells[i][j][k] == 1)
				{
					if (activecells[i-1][j][k] == 1) assert(std::copysign(1.0,sdfarray[i][j][k]) == std::copysign(1.0,sdfarray[i-1][j][k]));
					if (activecells[i+1][j][k] == 1) assert(std::copysign(1.0,sdfarray[i][j][k]) == std::copysign(1.0,sdfarray[i+1][j][k]));
					if (activecells[i][j-1][k] == 1) assert(std::copysign(1.0,sdfarray[i][j][k]) == std::copysign(1.0,sdfarray[i][j-1][k]));
					if (activecells[i][j+1][k] == 1) assert(std::copysign(1.0,sdfarray[i][j][k]) == std::copysign(1.0,sdfarray[i][j+1][k]));
					if (activecells[i][j][k-1] == 1) assert(std::copysign(1.0,sdfarray[i][j][k]) == std::copysign(1.0,sdfarray[i][j][k-1]));
					if (activecells[i][j][k+1] == 1) assert(std::copysign(1.0,sdfarray[i][j][k]) == std::copysign(1.0,sdfarray[i][j][k+1]));
				}
			}
		}
	}
}


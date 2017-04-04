#include "APSS.hpp"
#include "PS_storage.hpp"
#include "weightfunction.hpp"
#include "Nvector.hpp"
#include "fast_sweeping_method.hpp"
#include <fstream>
#include <algorithm>
#include <iostream>
#include <cassert>
#include <memory>
#include <vector>
#include <string>
#include <cstdlib>




void output_SDF_grid (const floatarray3D& sdfarray, unsigned int Nx, unsigned int Ny, unsigned int Nz, float x0, float y0, float z0, float dx, float dy, float dz, std::string filename);




int main()
{

	// Settings for code

	unsigned int Nx = 200;
	unsigned int Ny = 200;
	unsigned int Nz = 200;
	float x0 = -10;
	float y0 = -10;
	float z0 = -10;
	float dx = 0.1;
	float dy = 0.1;
	float dz = 0.1;
	std::string filename = "orion_nofbc";
	int extent = 1;

	
	std::cout << "Starting SDF generator on file: " << filename << std::endl;


	std::shared_ptr<weight_fn_base> wf = std::make_shared<polynomial_weight>();
	std::shared_ptr<PS_storage_3D> PS = std::make_shared<PS_storage_3D>(Nx,Ny,Nz,x0,y0,z0,dx,dy,dz,wf);
	APSS kenneth (PS, true);


	PS->load_point_set(filename + ".PN");
	PS->output_point_set(filename + ".dat");



	// Set size and intialise 3D arrays

	floatarray3D sdfarray;		// See fast_sweeping_method.hpp for type definition
	intarray3D activecells;
	sdfarray.resize(Nx);
	activecells.resize(Nx);
	int UBOUND = 1000000000;
	assert(UBOUND > extent);
	float UNSETVAL = 1e20;
	for (unsigned int i=0; i<Nx; ++i)
	{
		sdfarray[i].resize(Ny);
		activecells[i].resize(Ny);
		for (unsigned int j=0; j<Ny; ++j)
		{
			sdfarray[i][j].resize(Nz);
			activecells[i][j].resize(Nz);
			for (unsigned int k=0; k<Nz; ++k)
			{
				sdfarray[i][j][k] = UNSETVAL;
				activecells[i][j][k] = UBOUND;
			}
		}
	}
	
	

	// Decide on cells in which to compute SDF

	for (int i=0; i<int(Nx); ++i)
	{
		for (int j=0; j<int(Ny); ++j)
		{
			for (int k=0; k<int(Nz); ++k)
			{
				if (!(*PS).OParray[i][j][k].empty())
				{
					for (int a = std::max(0,i-extent); a<=std::min(int(Nx-1),i+extent); ++a)
					{
						for (int b = std::max(0,j-extent); b<=std::min(int(Ny-1),j+extent); ++b)
						{
							for (int c = std::max(0,k-extent); c<=std::min(int(Nz-1),k+extent); ++c)
							{
								activecells[a][b][c] = std::min(std::max({ abs(i-a), abs(j-b), abs(k-c)}), activecells[a][b][c]);
							}
						}
					}
				}
			}
		}

		std::cout << "\rSetting active cells.. " << (float(i)/Nx)*100 << "% done." << std::flush;
	}
	std::cout << "\rActive cells set.                                   " << std::endl;



	// Compute SDF in all cells within 'extent' of interface using APSS method
	
	float x, y, z, h1, h0 = 0.249999*std::min({dx, dy, dz});
	Nvector pos (3,0.0);

	for (unsigned int i=0; i<Nx; ++i)
	{
		x = x0 + 0.5*dx + i*dx;

		for (unsigned int j=0; j<Ny; ++j)
		{
			y = y0 + 0.5*dy + j*dy;

			for (unsigned int k=0; k<Nz; ++k)
			{
				if (activecells[i][j][k] < UBOUND)
				{

					z = z0 + 0.5*dz + k*dz;
					h = (activecells[i][j][k]*8.0 + 1.0)*h0;

					pos.data[0] = x;
					pos.data[1] = y;
					pos.data[2] = z;

					sdfarray[i][j][k] = kenneth.evaluate_surface(pos, h);
				}
			}
		}
		
		std::cout << "\rComputing SDF in active cells.. " << (float(i)/Nx)*100 << "% done." << std::flush;
	}
	std::cout << "\rAPSS computation complete.                                   " << std::endl;



	// Single iteration of the fast sweeping method should now by sufficient to set signed distance values in the rest of the domain

	std::cout << "Starting marching procedure." << std::endl;

	fast_sweeping_method (sdfarray, activecells, Nx, Ny, Nz, dx, dy, dz, UNSETVAL, UBOUND);
	
	std::cout << "\rAll computation complete.                                                         " << std::endl;





	output_SDF_grid(sdfarray,Nx,Ny,Nz,x0,y0,z0,dx,dy,dz,filename);

	std::cout << "Code complete!" << std::endl;

	return 0;
}


void output_SDF_grid (const floatarray3D& sdfarray, unsigned int Nx, unsigned int Ny, unsigned int Nz, float x0, float y0, float z0, float dx, float dy, float dz, std::string filename)
{
	/*
	 *	Store the 3D SDF in an ASCII VTK file for viewing in Visit
	 */

	std::ofstream vtkfile;
	vtkfile.open(filename + ".vtk");

	vtkfile << "# vtk DataFile Version 3.0\n";
	vtkfile << "vtkfile\n";
	vtkfile << "ASCII\n";
	vtkfile << "DATASET STRUCTURED_POINTS\n";
	vtkfile << "DIMENSIONS " <<Nx<< " " <<Ny<< " " <<Nz<< "\n";
	vtkfile << "ORIGIN " <<x0<<" "<<y0<<" " <<z0<<"\n";
	vtkfile << "SPACING "<<dx<<" "<<dy<<" "<<dz<<"\n";
	vtkfile << "POINT_DATA " << (Nx*Ny*Nz) << "\n";
	vtkfile << "SCALARS scalars float 1\n";
	vtkfile << "LOOKUP_TABLE default\n";

	for (unsigned int z=0; z<Nz; z++)
	{  
		for (unsigned int y=0; y<Ny; y++)
		{  
			for (unsigned int x=0; x<Nx; x++)
			{   
				vtkfile << sdfarray[x][y][z] << " ";
			}
		vtkfile << "\n";
		}
	vtkfile << "\n";
	std::cout << "\rWriting SDF to file.. " << (float(z)/Nz)*100 << "% done." << std::flush;
	}
	
	std::cout << "\rVTK output complete.                                   " << std::endl;
	vtkfile.close();
}

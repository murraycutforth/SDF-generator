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




void output_SDF_grid (const floatarray3D& sdfarray, const intarray3D& activecells, unsigned int Nx, unsigned int Ny, unsigned int Nz, float x0, float y0, float z0, float dx, float dy, float dz, std::string filename);




int main()
{

	// Settings for code

	unsigned int Nx = 300;
	unsigned int Ny = 300;
	unsigned int Nz = 300;
	float x0 = -10;
	float y0 = -10;
	float z0 = -10;
	float dx = 0.06666;
	float dy = 0.06666;
	float dz = 0.06666;
	std::string filename = "/local/data2/public/mcc74/APSS_STL/orion_nofbc";
	int extent = 1;

	
	std::cout << "Starting SDF generator on file: " << filename << std::endl;
	std::cout << "~" << std::endl;


	std::shared_ptr<weight_fn_base> wf = std::make_shared<polynomial_weight>();
	std::shared_ptr<PS_storage_3D> PS = std::make_shared<PS_storage_3D>(Nx,Ny,Nz,x0,y0,z0,dx,dy,dz,wf);
	APSS kenneth (PS, true);


	PS->load_point_set(filename + ".PN");
	std::cout << "~" << std::endl;



	// Set size and intialise 3D arrays

	floatarray3D sdfarray;		// See fast_sweeping_method.hpp for type definition
	intarray3D activecells;
	sdfarray.resize(Nx);
	activecells.resize(Nx);
	int UBOUND = 2;
	assert(UBOUND > extent);
	float UNSETVAL = 1e10;
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
	
	

	// Tag cells in narrow band around cells which contain marker particles
	
	int numAC0 = 0;

	for (int i=0; i<int(Nx); ++i)
	{
		for (int j=0; j<int(Ny); ++j)
		{
			for (int k=0; k<int(Nz); ++k)
			{
				if (!(*PS).OParray[i][j][k].empty())
				{
					numAC0++;
					
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
	std::cout << "Using this mesh there are " << numAC0 << " cells which contain marker particles." << std::endl;
	std::cout << "~" << std::endl;



	// Compute SDF in all tagged cells using APSS method
	
	float x, y, z, h, h0 = 0.249999*std::min({dx, dy, dz});
	Nvector pos (3,0.0);

	for (unsigned int i=0; i<Nx; ++i)
	{
		x = x0 + 0.5*dx + i*dx;

		for (unsigned int j=0; j<Ny; ++j)
		{
			y = y0 + 0.5*dy + j*dy;

			for (unsigned int k=0; k<Nz; ++k)
			{
				if (activecells[i][j][k] != UBOUND)
				{

					z = z0 + 0.5*dz + k*dz;
					h = (activecells[i][j][k]*8.0 + 1.0)*h0;

					pos.data[0] = x;
					pos.data[1] = y;
					pos.data[2] = z;

					sdfarray[i][j][k] = kenneth.evaluate_surface(pos, h);
					assert(sdfarray[i][j][k] != UNSETVAL);
				}
			}
		}
		
		std::cout << "\rComputing SDF in active cells.. " << (float(i)/Nx)*100 << "% done." << std::flush;
	}
	std::cout << "\rAPSS computation complete.                                   " << std::endl;
	std::cout << "~" << std::endl;
	
	
	
	// Prepare the activecells array to close the zero level set surface
	
	assert(extent == 1); // Need to do nothing for now if this is true
	
	
	// Single fast sweeping iteration should be sufficient to guarantee closed surface?
	
	std::cout << "Starting fast sweeping procedure to close the zero level set surface." << std::endl;
	
	fast_sweeping_method_closesurface (sdfarray, activecells, kenneth, Nx, Ny, Nz, dx, dy, dz, x0, y0, z0, UBOUND, UNSETVAL);
	
	std::cout << "\rSurface closed.                                 " << std::endl << "~" << std::endl;


	// Single iteration of the fast sweeping method should now by sufficient to set signed distance values in the rest of the domain

	std::cout << "Starting fast sweeping procedure to set value in bulk cells." << std::endl;

	fast_sweeping_method_fillbulkcells (sdfarray, activecells, Nx, Ny, Nz, dx, dy, dz, UNSETVAL, UBOUND);
	
	std::cout << "Sweeps complete." << std::endl << "~" << std::endl;
	
	std::cout << "\rAll computation complete.                                                         " << std::endl;





	output_SDF_grid(sdfarray,activecells,Nx,Ny,Nz,x0,y0,z0,dx,dy,dz,filename);

	std::cout << "Code complete!" << std::endl;

	return 0;
}


void output_SDF_grid (const floatarray3D& sdfarray, const intarray3D& activecells, unsigned int Nx, unsigned int Ny, unsigned int Nz, float x0, float y0, float z0, float dx, float dy, float dz, std::string filename)
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
	vtkfile << "SCALARS sdf float 1\n";
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

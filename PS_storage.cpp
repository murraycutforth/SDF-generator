#include "PS_storage.hpp"
#include "Nvector.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <cassert>
#include <algorithm>
#include <limits>


weighted_point :: weighted_point (float w, Nvector p, Nvector n)
:
	w (w),
	p (p),
	n (n)
{}


oriented_point :: oriented_point (Nvector p, Nvector n)
:
	p (p),
	n (n)
{}


int cellindex_x (float x0, float dx, float x)
{
	return static_cast<int>((x - x0)/dx);
}


int cellindex_y (float y0, float dy, float y)
{
	return static_cast<int>((y - y0)/dy);
}


int cellindex_z (float z0, float dz, float z)
{
	return static_cast<int>((z - z0)/dz);
}




PS_storage_3D :: PS_storage_3D (
	
	unsigned int Nx, 
	unsigned int Ny, 
	unsigned int Nz,
	float x0,
	float y0,
	float z0,
	float dx,
	float dy,
	float dz,
	std::shared_ptr<weight_fn_base> wf
)
:
	Nx	(Nx),
	Ny	(Ny),
	Nz	(Nz),
	x0	(x0),
	y0	(y0),
	z0	(z0),
	dx	(dx),
	dy	(dy),
	dz	(dz),
	OParray	(),
	wf	(wf)
{
	OParray.resize(Nx);
	for (unsigned int i=0; i<Nx; ++i)
	{
		OParray[i].resize(Ny);
		for (unsigned int j=0; j<Ny; ++j)
		{
			OParray[i][j].resize(Nz);
		}
	}
}



void PS_storage_3D :: load_point_set (std::string filename)
{
	std::ifstream infile (filename);
	
	// Assume that every line consists of 6 numbers
	
	float nx, ny, nz, x, y, z;
	Nvector normal (3, 0.0);
	Nvector position (3, 0.0);
	std::vector<oriented_point> temp_point_set;
	
	while (infile >> nx >> ny >> nz >> x >> y >> z)
	{
		normal.data[0] = nx;
		normal.data[1] = ny;
		normal.data[2] = nz;
		position.data[0] = x;
		position.data[1] = y;
		position.data[2] = z;


		// Partition the point cloud according to grid cells

		int i = cellindex_x(x0,dx,x);
		int j = cellindex_y(y0,dy,y);
		int k = cellindex_z(z0,dz,z);

		if (i >= 0 && j >= 0 && k >= 0
			&& i < int(Nx) && j < int(Ny) && k < int(Nz))
		{
			OParray[i][j][k].push_back(oriented_point(position, normal));
		}
	}

	infile.close();

	std::cout << "Point set loaded from file." << std::endl;
}



void PS_storage_3D :: get_weighted_points (

	const Nvector& x, 
	std::vector<weighted_point>& wps, 
	const float h
) const
{
	/*
	 * In current implementation, we search through all points
	 * in local grid cells, and calculate the weighting of each one.
	 */
	
	
	wps.clear();
	float wi;

	unsigned int i_min, i_max, j_min, j_max, k_min, k_max;

	int imint, imaxt, jmint, jmaxt, kmint, kmaxt;

	imint = cellindex_x(x0,dx,x.data[0]-h);
	imaxt = cellindex_x(x0,dx,x.data[0]+h);
	jmint = cellindex_y(y0,dy,x.data[1]-h);
	jmaxt = cellindex_y(y0,dy,x.data[1]+h);
	kmint = cellindex_z(z0,dz,x.data[2]-h);
	kmaxt = cellindex_z(z0,dz,x.data[2]+h);

	i_min = std::max(0,imint);
	i_max = std::min(int(Nx)-1,imaxt);
	j_min = std::max(0,jmint);
	j_max = std::min(int(Ny)-1,jmaxt);
	k_min = std::max(0,kmint);
	k_max = std::min(int(Nz)-1,kmaxt);

	for (unsigned int i=i_min; i<=i_max; ++i)
	{
		for (unsigned int j=j_min; j<=j_max; ++j)
		{
			for (unsigned int k=k_min; k<=k_max; ++k)
			{
				for (const oriented_point& OP : OParray[i][j][k])
				{
					wi = wf->weight((OP.p - x).magnitudesq() / (h*h));

					if (wi > 10.0*std::numeric_limits<float>::epsilon())
					{
						wps.push_back(weighted_point(wi, OP.p, OP.n));
					}
				}
			}
		}
	}
}



void PS_storage_3D :: output_point_set (std::string filename)
{
	std::ofstream outfile;
	outfile.open(filename);
	
	for (unsigned int i=0; i<Nx; ++i)
	{
		for (unsigned int j=0; j<Ny; ++j)
		{
			for (unsigned int k=0; k<Nz; ++k)
			{
				for (const oriented_point& OP : OParray[i][j][k])
				{
					for (float coord : OP.p.data)
					{
						outfile << coord << " ";
					}
					outfile << std::endl;
				}
			}
		}
	}
	
	std::cout << "\r" << "Point set output complete" << std::endl;
	outfile.close();
}






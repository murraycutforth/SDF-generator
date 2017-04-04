/*
 * DESCRIPTION
 * 	This class is used to store a set of points with normals. The get_all_weighted_points function
 *	is called by the APSS function in order to find the weights of all points in the cloud with respect
 *	to some evaluation position. 
 * 
 * AUTHOR
 * 	Murray Cutforth
 * 
 * Date
 * 	20/11/2015
 */


#ifndef PSSTORAGE_H
#define PSSTORAGE_H


#include "Nvector.hpp"
#include "weightfunction.hpp"
#include <vector>
#include <string>
#include <memory>




struct weighted_point {
	
	float w;		// Weighting of this point
	Nvector p;		// Position of point
	Nvector n;		// Unit normal vector of point
	
	weighted_point (float w, Nvector p, Nvector n);
};


struct oriented_point {
	
	Nvector p;		// Position of point
	Nvector n;		// Unit normal vector of point
	
	oriented_point (Nvector p, Nvector n);
};


typedef std::vector<std::vector<std::vector<std::vector<oriented_point> > > > OParray4D;


int cellindex_x (float x0, float dx, float x);

int cellindex_y (float y0, float dy, float y);

int cellindex_z (float z0, float dz, float z);




class PS_storage_3D {

	public:
		
	unsigned int Nx; 
	unsigned int Ny; 
	unsigned int Nz;
	float x0;
	float y0;
	float z0;
	float dx;
	float dy;
	float dz;

	OParray4D OParray;

	std::shared_ptr<weight_fn_base> wf;

	PS_storage_3D (	
	
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
	);


	void load_point_set (std::string filename);

	void get_weighted_points (const Nvector& x, std::vector<weighted_point>& wps, const float h) const;

	void output_point_set (std::string filename);
};























class PS_storage_base {
	
	public:
		
	std::shared_ptr<weight_fn_base> wf;				// Pointer to weight function object
	
	
	// Constructors
	
	PS_storage_base ();
	
	PS_storage_base (std::shared_ptr<weight_fn_base> wf);
	
	
	// Function which initialises the point_set with data from file
	
	virtual void fill_point_set (std::string filename) =0;
	
	
	// Return all points with nonzero weighting from x given the length scale h

	virtual void get_all_weighted_points (const Nvector& x, std::vector<weighted_point>& wps, const float h) const =0;


	// Test out some predefined point clouds

	virtual void test_fill_point_set (int option, int N) =0;
	
	
	// Print out all marker positions
	
	virtual void output_markers (std::string filename) =0;
};




class PS_storage_naive : public PS_storage_base {

	/*
	 *	This class stores all points in the cloud in one disordered vector. Every time we want to get the weighting of each point
	 *	with respect to some evaluation point, this class iterates over every single point in the point cloud and computes their weights.
	 */
	
	public:

	std::vector<oriented_point> point_set;
	
	
	PS_storage_naive ();
		
	PS_storage_naive (std::shared_ptr<weight_fn_base> wf);
	
	
	void fill_point_set (std::string filename);
	
	void test_fill_point_set (int option, int N);
	
	void get_all_weighted_points (const Nvector& x, std::vector<weighted_point>& wps, const float hx) const;	
	
	void output_markers (std::string filename);
	
};




#endif

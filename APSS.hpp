/*
 * DESCRIPTION
 * 	This class implements the algebraic point set surface method to
 * 	compute a signed distance function from a set of oriented points
 * 	which lie on a surface. The method uses the class in "Nvector.hpp"
 * 	to work independently of dimension.
 * 
 * AUTHOR
 * 	Murray Cutforth
 * 
 * Date
 * 	4/3/2016
 */


#ifndef APSSMAIN_H
#define APSSMAIN_H


#include "PS_storage.hpp"
#include "Nvector.hpp"
#include <memory>	


class APSS {
	
	public:
	
	std::shared_ptr<PS_storage_3D> PS_ptr;			// Storage of point set cloud
	
	bool negativeinterior;					// If this is true, we use convention that the interior of the surface has negative distance function
	
	APSS (std::shared_ptr<PS_storage_3D> ps_ptr, bool negativeinterior);
	
	
	// Compute the APSS at a given point. Return the value of the APSS (apssval), and the signed distance (sdval).
	
	float evaluate_surface (const Nvector& x, float h0) const;
};

#endif

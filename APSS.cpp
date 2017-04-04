#include "APSS.hpp"
#include <vector>
#include <limits>
#include <iostream>
#include <cassert>
#include <cmath>


APSS :: APSS (std::shared_ptr<PS_storage_3D> ps_ptr, bool negativeinterior)
:
	PS_ptr (ps_ptr),
	negativeinterior (negativeinterior)
{}


float APSS :: evaluate_surface (const Nvector& x, float h0) const
{
	
	// First extract all points with nonzero weighting and store in wps
	
	const float beta = 1.0;
	std::vector<weighted_point> wps;
	double h = h0;
	
		
	// Require at least 1 marker - double length scale until we find nonzero weighted points
	
	while (wps.empty())
	{
		PS_ptr->get_weighted_points(x,wps,h);
		h *= 2.0;
	}

//	PS_ptr->get_weighted_points(x,wps,h);
//	if (wps.empty())
//	{
//		std::cout << std::endl << "Empty point set! Using h = " << h << std::endl;
//		std::cout << "pos = " << x.data[0] << "," << x.data[1] << "," << x.data[2] << std::endl;
//	}
//	assert(!wps.empty());
	
	
	// Initialise storage for various sums
	
	float Sw = 0;
	float Swnp = 0;
	float Swpp = 0;
	float SwnSwp = 0;
	float SwpSwp = 0;
	Nvector Swp_ (x.get_dimension(), 0.0);
	Nvector Swn_ (x.get_dimension(), 0.0);
	
	
	// Now iterate over all points with nonzero weighting to compute necessary quantities
	
	for (const weighted_point& wpi : wps)
	{
		Sw += wpi.w;
		
		Swp_ += wpi.w * wpi.p;
	}
	
	for (const weighted_point& wpi : wps)
	{
		
		Swnp += wpi.w * wpi.n.inner_product(wpi.p);
		
		Swpp += wpi.w * wpi.p.inner_product(wpi.p);
		
		SwnSwp += wpi.w * wpi.n.inner_product(Swp_);
			
		SwpSwp += wpi.w * wpi.p.inner_product(Swp_);
	
		Swn_ += wpi.w * wpi.n;
		
	}
	
	assert(Sw != 0.0);
	
	
	// Calculate the coefficients of u - accounting for degenerate case where sphere is actually a plane
	
	float u_d1;
	
	if (fabs(Sw * Swpp - SwpSwp) < 10.0*std::numeric_limits<float>::epsilon())
	{
		u_d1 = 0.0;
	}
	else
	{
		u_d1 = beta * 0.5 * ((Sw * Swnp - SwnSwp) / (Sw * Swpp - SwpSwp));
	}
		
	Nvector u_mid ((1.0/Sw) * (Swn_ - (2.0 * u_d1 * Swp_)));
	
	float u0 = - (1.0/Sw) * (u_mid.inner_product(Swp_) + u_d1 * Swpp);
	
	
	// Estimate the distance from x to the surface of the fitted sphere

	float sdval;
		
	if (fabs(u_d1) > 0.0)
	{
	
		Nvector c (-0.5*(1.0/u_d1)*u_mid);
	
		float r = sqrt(c.inner_product(c) - (u0/u_d1));
		
		sdval = fabs((c - x).magnitude() - r);
	}
	else
	{
		
		// In this case the sphere has degenerated to a plane
		
		sdval = fabs(x.inner_product(u_mid) + u0) / u_mid.magnitude();
	}


	// Compute value of algebraic surface at this point

	float apssval;

	apssval = u0 + x.inner_product(u_mid) + x.inner_product(x)*u_d1;

	if (!negativeinterior) apssval *= -1.0;

	
	// Give the appropriate sign to the distance
	
	sdval = std::copysign(sdval, apssval);

	return sdval;
}


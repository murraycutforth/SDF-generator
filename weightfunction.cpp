#include "weightfunction.hpp"
#include <cassert>
#include <cmath>


float polynomial_weight :: weight (const float x) const
{
	assert(x >= 0.0);
	
	if (x < 1.0)
	{
		return (1.0 - x)*(1.0 - x)*(1.0 - x)*(1.0 - x);
	}
	else
	{
		return 0.0;
	}
}



float narrow_polynomial_weight :: weight (const float x) const
{
	assert(x >= 0.0);
	
	if (x < 1.0)
	{
		return std::pow((1.0 - x), 8);
	}
	else
	{
		return 0.0;
	}
}

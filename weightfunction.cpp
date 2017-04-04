#include "weightfunction.hpp"
#include <cassert>


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

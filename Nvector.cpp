#include "Nvector.hpp"
#include <cmath>
#include <cassert>
#include <algorithm>


Nvector :: Nvector ()
:
	data ()
{}


Nvector :: Nvector (unsigned int n, float val)
:
	data (n,val)
{}


Nvector :: Nvector (const Nvector& other)
:
	data (other.data)
{}



unsigned int Nvector :: get_dimension () const
{
	return data.size();
}


float Nvector :: magnitude () const
{
	float ss = 0.0;
	
	for (float x : data)
	{
		ss += x*x;
	}
	
	return sqrt(ss);
}


float Nvector :: magnitudesq () const
{
	float ss = 0.0;
	
	for (float x : data)
	{
		ss += x*x;
	}
	
	return ss;
}


float Nvector :: inner_product (const Nvector& other) const
{
	assert(get_dimension() == other.get_dimension());
	
	float sum = 0.0;
	
	for (unsigned int i=0; i<get_dimension(); i++)
	{
		sum += data[i]*other.data[i];
	}
	
	return sum;
}




void swap (Nvector& first, Nvector& second)
{
	std::swap(first.data,second.data);
}


Nvector& Nvector :: operator= (Nvector other)
{
	swap(*this,other);
	
	return *this;
}


Nvector& Nvector :: operator+= (const Nvector& rhs)
{
	assert(this->get_dimension() == rhs.get_dimension());
	
	for (unsigned int i=0; i<this->get_dimension(); i++)
	{
		data[i] += rhs.data[i];
	}
	
	return *this;
}	
	
	
Nvector& Nvector :: operator-= (const Nvector& rhs)
{
	assert(this->get_dimension() == rhs.get_dimension());
	
	for (unsigned int i=0; i<this->get_dimension(); i++)
	{
		data[i] -= rhs.data[i];
	}
	
	return *this;
}


Nvector& Nvector :: operator*= (float rhs)
{
	for (float& x : data)
	{
		x *= rhs;
	}
	
	return *this;
}


Nvector operator+ (Nvector lhs, const Nvector& rhs)
{
	assert(lhs.get_dimension() == rhs.get_dimension());
	
	lhs += rhs;
	return lhs;
}


Nvector operator- (Nvector lhs, const Nvector& rhs)
{
	assert(lhs.get_dimension() == rhs.get_dimension());
	
	lhs -= rhs;
	return lhs;
}


Nvector operator* (Nvector lhs, float rhs)
{
	lhs *= rhs;
	return lhs;
}


Nvector operator* (float lhs, Nvector rhs)
{
	rhs *= lhs;
	return rhs;
}

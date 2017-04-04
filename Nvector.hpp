/*
 * DESCRIPTION
 * 	This class describes an n-dimensional vector, with overloads for
 *	a few of the necessary mathematical operations which are used by the
 *	APSS method.
 * 
 * AUTHOR
 * 	Murray Cutforth
 * 
 * Date
 * 	4/3/2016
 */


#ifndef NVECTOR_H
#define NVECTOR_H


#include <vector>


class Nvector {
	
	public:
	
	std::vector<float> data;				// Storage for components of vector
	
	
	// Constructors

	Nvector ();
	
	Nvector (unsigned int n, float val);
	
	Nvector (const Nvector& other);
	

	// Mathematic functions
	
	unsigned int get_dimension () const;
	
	float magnitude () const;
	
	float magnitudesq () const;
	
	float inner_product (const Nvector& other) const;
	
	
	Nvector& operator= (Nvector other);
	
	Nvector& operator+= (const Nvector& rhs);
	
	Nvector& operator-= (const Nvector& rhs);
	
	Nvector& operator*= (float rhs);

		
};


void swap (Nvector& first, Nvector& second);

Nvector operator+ (Nvector lhs, const Nvector& rhs);

Nvector operator- (Nvector lhs, const Nvector& rhs);

Nvector operator* (Nvector lhs, float rhs);

Nvector operator* (float lhs, Nvector rhs);


#endif

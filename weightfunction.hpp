/*
 * DESCRIPTION
 * 	These classes define the weight function used in the APSS which is used to weight each point in the surface fitting.
 * 
 * AUTHOR
 * 	Murray Cutforth
 * 
 * Date
 * 	23/11/2015
 */


#ifndef APSSWEIGHTFUNCTION_H
#define APSSWEIGHTFUNCTION_H


class weight_fn_base {
	
	public:
	
	virtual float weight (float x) const =0;
};


class polynomial_weight : public weight_fn_base {
	
	public:
	
	float weight (float x) const;
};

#endif

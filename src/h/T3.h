/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     T3 Element Extension                                                  */
/*****************************************************************************/

#pragma once

#include "Element.h"

using namespace std;

//! 3-node triangular element class (Constant Strain Triangle, CST)
class CT3 : public CElement
{
public:

//!	Constructor
	CT3();

//!	Desconstructor
	~CT3();

//!	Read element data from stream Input
	virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);

//!	Write element data to stream
	virtual void Write(COutputter& output);

//!	Generate location matrix
	virtual void GenerateLocationMatrix();

//!	Calculate element stiffness matrix
	virtual void ElementStiffness(double* Matrix);

//!	Calculate element stress
	virtual void ElementStress(double* stress, double* Displacement);

//!	Return the size of the element stiffness matrix (upper triangular)
	virtual unsigned int SizeOfStiffnessMatrix() { return 21; }
};
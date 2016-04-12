#ifndef _ILUBase_
#define _ILUBase_

#include <iostream>
#include "SparseMatrix.h"

class ILUBase
{
public:
	ILUBase(int p = 0) {};
	ILUBase(SparseMatrixCRS &A, SparseMatrixCRS &M, int p) {} ;
	virtual void Compute(SparseMatrixCRS &A, SparseMatrixCRS &M, int p) = 0;
	virtual void Compute() = 0;
	SparseMatrixCRS *GetA() { return mA; }
	SparseMatrixCRS *GetM() { return mM; }
	virtual ~ILUBase() {} ;
protected:
	SparseMatrixCRS *mA;
	SparseMatrixCRS *mM;
	int mP;
};

#endif
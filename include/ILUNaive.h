#ifndef _ILUNaive_
#define _ILUNaive_


#include <iostream>
#include "SparseMatrix.h"
#include "ILUBase.h"

class ILUNaive : public ILUBase
{
public:
	ILUNaive(int p = 0)	{ mP = p;}
	ILUNaive(SparseMatrixCRS &A, SparseMatrixCRS &M, int p) {} ;
	bool CheckAinM(SparseMatrixCRS &A, SparseMatrixCRS &M);
	bool ILUNaive::CheckInverse(SparseMatrixCRS &A, SparseMatrixCRS &M);
	virtual void Compute(SparseMatrixCRS &A, SparseMatrixCRS &M, int p);
	virtual void Compute() {};
	SparseMatrixCRS *GetA() { return mA; }
	SparseMatrixCRS *GetM() { return mM; }
	virtual ~ILUNaive() {};
private:
	SparseMatrixCRS *mA;
	SparseMatrixCRS *mM;
	int mP;
};

#endif
#ifndef _ILUOpenMP_
#define _ILUOpenMP_


#include <iostream>
#include "SparseMatrix.h"

class ILU
{
public:
	ILU(int p = 0)	{ mP = p;}
	ILU(SparseMatrixCRS &A, SparseMatrixCRS &M, int p) {} ;
	bool CheckAinM(SparseMatrixCRS &A, SparseMatrixCRS &M);
	bool CheckInverse(SparseMatrixCRS &A, SparseMatrixCRS &M);
	bool isCorrectMatrix(SparseMatrixCRS &A);
	void Compute(SparseMatrixCRS &A, SparseMatrixCRS &M, int p, int flags);
	void Compute() {};
	SparseMatrixCRS *GetA() { return mA; }
	SparseMatrixCRS *GetM() { return mM; }
	~ILU() {};
private:
	SparseMatrixCRS *mA;
	SparseMatrixCRS *mM;
	int mP;
};

#endif
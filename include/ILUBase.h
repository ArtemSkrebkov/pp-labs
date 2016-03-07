#include <iostream>
#include "SparseMatrix.h"

class ILUBase
{
public:
	ILUBase(int p);
	ILUBase(SparseMatrixBase *A, SparseMatrixBase *M, int p);
	virtual void Compute(SparseMatrixBase *A, SparseMatrixBase *M, int p) = 0;
	virtual void Compute() = 0;
	SparseMatrixBase *GetA() { return mA; }
	SparseMatrixBase *GetM() { return mM; }
	virtual ~ILUBase();
protected:
	SparseMatrixBase *mA;
	SparseMatrixBase *mM;
	int p;
};
#include <iostream>
#include "SparseMatrix.h"
#include "ILUBase.h"

class ILUTest
{
public:
	ILUTest(int p);
	ILUTest(SparseMatrixBase *A, SparseMatrixBase *M, int p);
	virtual void Compute(SparseMatrixBase *A, SparseMatrixBase *M, int p);
	virtual void Compute();
	virtual ~ILUTest();
};
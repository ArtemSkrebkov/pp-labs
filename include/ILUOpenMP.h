#ifndef _ILUOpenMP_
#define _ILUOpenMP_


#include <iostream>
#include "SparseMatrix.h"
#include "ILUBase.h"

class ILUOpenMP : public ILUBase
{
public:
	ILUOpenMP(int p = 0)	{ mP = p;}
	ILUOpenMP(SparseMatrixCRS &A, SparseMatrixCRS &M, int p) {} ;
	bool CheckAinM(SparseMatrixCRS &A, SparseMatrixCRS &M);
	bool CheckInverse(SparseMatrixCRS &A, SparseMatrixCRS &M);
	bool isCorrectMatrix(SparseMatrixCRS &A);
	virtual void Compute(SparseMatrixCRS &A, SparseMatrixCRS &M, int p);
	virtual void Compute() {};
	SparseMatrixCRS *GetA() { return mA; }
	SparseMatrixCRS *GetM() { return mM; }
	virtual ~ILUOpenMP() {};
private:
	SparseMatrixCRS *mA;
	SparseMatrixCRS *mM;
	int mP;
};

class ILUTBB : public ILUBase
{
public:
	ILUTBB(int p = 0)	{ mP = p;}
	ILUTBB(SparseMatrixCRS &A, SparseMatrixCRS &M, int p) {} ;
	bool CheckAinM(SparseMatrixCRS &A, SparseMatrixCRS &M);
	bool CheckInverse(SparseMatrixCRS &A, SparseMatrixCRS &M);
	bool isCorrectMatrix(SparseMatrixCRS &A);
	virtual void Compute(SparseMatrixCRS &A, SparseMatrixCRS &M, int p);
	virtual void Compute() {};
	SparseMatrixCRS *GetA() { return mA; }
	SparseMatrixCRS *GetM() { return mM; }
	virtual ~ILUTBB() {};
private:
	SparseMatrixCRS *mA;
	SparseMatrixCRS *mM;
	int mP;
};

class ILUCilk : public ILUBase
{
public:
	ILUCilk(int p = 0)	{ mP = p;}
	ILUCilk(SparseMatrixCRS &A, SparseMatrixCRS &M, int p) {} ;
	bool CheckAinM(SparseMatrixCRS &A, SparseMatrixCRS &M);
	bool CheckInverse(SparseMatrixCRS &A, SparseMatrixCRS &M);
	bool isCorrectMatrix(SparseMatrixCRS &A);
	virtual void Compute(SparseMatrixCRS &A, SparseMatrixCRS &M, int p);
	virtual void Compute() {};
	SparseMatrixCRS *GetA() { return mA; }
	SparseMatrixCRS *GetM() { return mM; }
	virtual ~ILUCilk() {};
private:
	SparseMatrixCRS *mA;
	SparseMatrixCRS *mM;
	int mP;
};


#endif
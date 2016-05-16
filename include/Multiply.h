#ifndef _MULTIPLY_
#define _MULTIPLY_

#include "Matrix.h"
#include "SparseMatrix.h"

void MultiplyFullMatrix(Matrix &a, Matrix &b, Matrix &c);

void MultiplyNaive(SparseMatrixCRS &A, SparseMatrixCRS &B, SparseMatrixCRS &C);
void MultiplyOpenMP(SparseMatrixCRS &A, SparseMatrixCRS &B, SparseMatrixCRS &C);
void MultiplyTBB(SparseMatrixCRS &A, SparseMatrixCRS &B, SparseMatrixCRS &C);
void MultiplyCilk(SparseMatrixCRS &A, SparseMatrixCRS &B, SparseMatrixCRS &C);

#endif
#include <gtest/gtest.h>
#include <iostream>
#include "SparseMatrix.h"
#include "ILUNaive.h"

using namespace std;

TEST(ILU, correct_ilu)
{
    // TO DO
	ILUNaive ilu(0);

	SparseMatrixCRS A("testdata\\A.txt");
	SparseMatrixCRS M;
	ilu.Compute(A, M, 0);

	SparseMatrixCRS cA("testdata\\A.txt");
	EXPECT_EQ(true, ilu.CheckAinM(cA, M));
	EXPECT_EQ(true, ilu.CheckInverse(cA, M));
}

TEST(ILU, GershgorinConditionNumber)
{
	SparseMatrixCRS A("testdata\\A.txt");
    EXPECT_EQ(3, A.GershgorinConditionNumber());
}


TEST(MATRIX, transpose)
{
	SparseMatrixCRS A("testdata\\AT.txt");
	SparseMatrixCRS AT = A.Transpose();
	SparseMatrixCRS ATT = AT.Transpose();
    EXPECT_EQ(A, ATT);
}

TEST(MATRIX, multiply)
{
	SparseMatrixCRS A("testdata\\A.txt");
	SparseMatrixCRS B("testdata\\A.txt");
	SparseMatrixCRS BT = B.Transpose();
	SparseMatrixCRS rightAnswer("testdata\\AA.txt");
	SparseMatrixCRS C;
	A.Multiply(A, B, C);
    EXPECT_EQ(rightAnswer, C);
}
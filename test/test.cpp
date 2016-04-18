#include <gtest/gtest.h>
#include <iostream>
#include "SparseMatrix.h"
#include "ILUNaive.h"

using namespace std;

TEST(ILU_PERF, simple)
{
    // TO DO
	ILUNaive ilu(0);

	SparseMatrixCRS A("C:\\Users\\Artem\\Documents\\PP_labs\\testdata\\cdde6.mtx");
	SparseMatrixCRS M;

	const clock_t t0 = clock(); // or gettimeofday or whatever
	ilu.Compute(A, M, 0);
	const clock_t t1 = clock();
	const double elapsedSec = (t1 - t0) / (double)CLOCKS_PER_SEC;

	cout << "elapsedSec = " << elapsedSec << endl;
	
	EXPECT_EQ(true, ilu.CheckAinM(A, M));
	EXPECT_EQ(true, ilu.CheckInverse(A, M));
}

TEST(ILU, correct_ilu)
{
    // TO DO
	ILUNaive ilu(0);

	SparseMatrixCRS A("C:\\Users\\Artem\\Documents\\PP_labs\\testdata\\A.txt");
	SparseMatrixCRS M;
	ilu.Compute(A, M, 0);

	SparseMatrixCRS cA("C:\\Users\\Artem\\Documents\\PP_labs\\testdata\\A.txt");
	EXPECT_EQ(true, ilu.CheckAinM(cA, M));
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

TEST(MATRIX, inverse)
{
	vector<vector<double>> M;
	M.resize(10);
	for (int i = 0;i<10; i++)
	{
		M[i].resize(10);
	}

}
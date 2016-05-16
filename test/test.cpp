#include <gtest/gtest.h>
#include <iostream>
#include "SparseMatrix.h"
#include "ILU.h"
#include "Multiply.h"

using namespace std;

TEST(ILU_PERF, simple)
{
	ILU ilu(0);
	SparseMatrixCRS A("/home/artem/workspace/build/pp-labs-build/bin/testdata/msc00726.mtx");

	SparseMatrixCRS M;
	ASSERT_EQ(true, ilu.isCorrectMatrix(A));
	const clock_t t0 = clock(); // or gettimeofday or whatever
	ilu.Compute(A, M, 0, 0);
	const clock_t t1 = clock();
	const double elapsedSec = (t1 - t0) / (double)CLOCKS_PER_SEC;

	cout << "elapsedSec = " << elapsedSec << endl;
	
	EXPECT_EQ(true, ilu.CheckAinM(A, M));
	EXPECT_EQ(true, ilu.CheckInverse(A, M));
}

TEST(ILU_PERF, openmp)
{
	ILU ilu(0);

	SparseMatrixCRS A("/home/artem/workspace/build/pp-labs-build/bin/testdata/msc00726.mtx");
	SparseMatrixCRS M;
	ASSERT_EQ(true, ilu.isCorrectMatrix(A));
	const clock_t t0 = clock(); // or gettimeofday or whatever
	ilu.Compute(A, M, 0, 1);
	const clock_t t1 = clock();
	const double elapsedSec = (t1 - t0) / (double)CLOCKS_PER_SEC;

	cout << "elapsedSec = " << elapsedSec << endl;
	
	EXPECT_EQ(true, ilu.CheckAinM(A, M));
	EXPECT_EQ(true, ilu.CheckInverse(A, M));
}

TEST(ILU_PERF, tbb)
{
	ILU ilu(0);

	SparseMatrixCRS A("/home/artem/workspace/build/pp-labs-build/bin/testdata/msc00726.mtx");
	SparseMatrixCRS M;
	ASSERT_EQ(true, ilu.isCorrectMatrix(A));
	const clock_t t0 = clock(); // or gettimeofday or whatever
	ilu.Compute(A, M, 0, 2);
	const clock_t t1 = clock();
	const double elapsedSec = (t1 - t0) / (double)CLOCKS_PER_SEC;

	cout << "elapsedSec = " << elapsedSec << endl;
	
	EXPECT_EQ(true, ilu.CheckAinM(A, M));
	EXPECT_EQ(true, ilu.CheckInverse(A, M));
}


TEST(ILU_PERF, cilk)
{
	ILU ilu(0);

	SparseMatrixCRS A("/home/artem/workspace/build/pp-labs-build/bin/testdata/msc00726.mtx");
	SparseMatrixCRS M;
	ASSERT_EQ(true, ilu.isCorrectMatrix(A));
	const clock_t t0 = clock(); // or gettimeofday or whatever
	ilu.Compute(A, M, 0, 3);
	const clock_t t1 = clock();
	const double elapsedSec = (t1 - t0) / (double)CLOCKS_PER_SEC;

	cout << "elapsedSec = " << elapsedSec << endl;
	
	EXPECT_EQ(true, ilu.CheckAinM(A, M));
	EXPECT_EQ(true, ilu.CheckInverse(A, M));
}

TEST(ILU, correct_ilu)
{
	ILU ilu(0);

	SparseMatrixCRS A("/home/artem/workspace/build/pp-labs-build/bin/testdata/A.txt");
	SparseMatrixCRS M;
	ilu.Compute(A, M, 0, 0);

	SparseMatrixCRS cA("/home/artem/workspace/build/pp-labs-build/bin/testdata/A.txt");
	EXPECT_EQ(true, ilu.CheckAinM(cA, M));
}

TEST(MATRIX, transpose)
{
	SparseMatrixCRS A("/home/artem/workspace/build/pp-labs-build/bin/testdata/AT.txt");
	SparseMatrixCRS AT = A.Transpose();
	SparseMatrixCRS ATT = AT.Transpose();
    EXPECT_EQ(A, ATT);
}

TEST(MATRIX, multiply)
{
	SparseMatrixCRS A("/home/artem/workspace/build/pp-labs-build/bin/testdata/A.txt");
	SparseMatrixCRS B("/home/artem/workspace/build/pp-labs-build/bin/testdata/A.txt");
	SparseMatrixCRS BT = B.Transpose();
	SparseMatrixCRS rightAnswer("/home/artem/workspace/build/pp-labs-build/bin/testdata/AA.txt");
	SparseMatrixCRS C;

	MultiplyNaive(A, B, C);
    EXPECT_EQ(rightAnswer, C);

	MultiplyOpenMP(A, B, C);
    EXPECT_EQ(rightAnswer, C);

	MultiplyTBB(A, B, C);
    EXPECT_EQ(rightAnswer, C);

    MultiplyCilk(A, B, C);
    EXPECT_EQ(rightAnswer, C);
}
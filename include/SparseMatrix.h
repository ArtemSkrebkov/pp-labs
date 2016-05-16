#ifndef _SPARSE_MATRIX_
#define _SPARSE_MATRIX_


#include <vector>
#include <fstream>
#include <cmath>

const double ZERO_IN_CRS = 0.00000001;

class SparseMatrixCRS
{
public:
    SparseMatrixCRS();
	SparseMatrixCRS(const SparseMatrixCRS &c);
    SparseMatrixCRS(const std::string &filename);
	SparseMatrixCRS Transpose();
	bool IsNonZero(size_t i, size_t j);
    double Get(int i, int j) const;
    void Set(int i, int j, double val) ;
    ~SparseMatrixCRS();
    void Print( size_t outSize = 0, size_t start_i = 0, size_t start_j = 0);
    void ReadFromMtx(const std::string filename);

	inline bool operator==(const SparseMatrixCRS& rhs)
	{
		bool result = true;
		if (this->mN != rhs.mN || this->mNZ != rhs.mNZ)
		{
			result = false;
		}
		else
		{
			for (size_t i = 0; i < rhs.mN && result; i++)
			{
				for (size_t j = 0; j < rhs.mN && result; j++)
				{
					if (fabs(this->Get(i, j) - rhs.Get(i, j)) > ZERO_IN_CRS)
					{
						result = false;
						break;
					}
				}
			}
		}

		return result;
	}
	inline bool operator!=(const SparseMatrixCRS& rhs){ return !(*this == rhs); }

	friend inline bool operator==(const SparseMatrixCRS& lhs, const SparseMatrixCRS& rhs);
public:
	void InitializeMatrix(size_t N, size_t NZ);
    size_t mN, mNZ;

    std::vector<double> mValues;
    std::vector<size_t> mCol;
    std::vector<size_t> mRowIndex;
};

inline bool operator==(const SparseMatrixCRS& lhs, const SparseMatrixCRS& rhs)
	{
		bool result = true;
		if (lhs.mN != rhs.mN || lhs.mNZ != rhs.mNZ)
		{
			result = false;
		}
		else
		{
			for (size_t i = 0; i < rhs.mN; i++)
			{
				for (size_t j = 0; j < rhs.mN; j++)
				{
					if (lhs.Get(i, j) != rhs.Get(i, j))
					{
						result = false;
						break;
					}
				}
			}
		}

		return result;
	}
	

#endif
#ifndef MATRICES_H
#define MATRICES_H

#include "engine/math/hcVec.h"

#ifdef GUI
#include <extern/assimp/scene.h>
#endif

#include <stdlib.h>

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          Matrix
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

/*! \brief implementation of the mathematical matrix construct, only float values supported so far
 *
 * 	Stores entries of a matrix in row-major ordering in array content
 */
template<uint rows, uint cols, class T>
class Matrix{
public:

	T content[rows * cols];												/*!< \brief contains the entries of the matrix in row-major ordering 	*/

#ifdef __NVCC__
__host__ __device__
#endif
	Matrix();                                 							/*!< \brief std constructor 					*/

template<class S>
#ifdef __NVCC__
__host__ __device__
#endif
	Matrix(const Matrix<rows, cols, S> &other);   						/*!< \brief copy constructor 					*/

#ifdef __NVCC__
__host__ __device__
#endif
	Matrix(const Vec<rows, T> &v1, const Vec<rows, T> &v2);           	/*!< \brief create column-matrix 	*/


#ifdef __NVCC__
__host__ __device__
#endif
	Matrix(const Vec<rows, T> &v1, const Vec<rows, T> &v2, const Vec<rows, T> &v3);    /*!< \brief create column-matrix 	*/

#ifdef __NVCC__
__host__ __device__
#endif
	~Matrix();                                 								/*!< \brief destructor 					*/

template<class S>
#ifdef __NVCC__
__host__ __device__
#endif
	Matrix &operator=(const Matrix<rows, cols, S> &other);

#ifdef __NVCC__
__host__ __device__
#endif
    Matrix &operator*=(T scale);
	/*!< \brief scales the matrix                                                           */

#ifdef __NVCC__
__host__ __device__
#endif
    Matrix &operator/=(T scale);
        /*!< \brief scales the matrix                                                           */

#ifdef __NVCC__
__host__ __device__
#endif
    Matrix &operator+=(const Matrix<rows, cols, T> &other);
        /*!< \brief adds another matrix to this                                                 */

#ifdef __NVCC__
__host__ __device__
#endif
    Matrix &operator-=(const Matrix<rows, cols, T> &other);
        /*!< \brief adds another matrix to this                                                 */

#ifdef __NVCC__
__host__ __device__
#endif
    Vec<rows, T> operator*(const Vec<cols, T> &vec);
        /*!< \brief this * vec                                                                  */

#ifdef __NVCC__
__host__ __device__
#endif
	T  &operator()(uint i, uint j);
		/*!< \brief access to element in matrix at row j, column j								*/

#ifdef __NVCC__
__host__ __device__
#endif
	T operator[](uint i) const;
		/*!< \brief access to content in row-major ordering										*/

#ifdef __NVCC__
__host__ __device__
#endif
	T getElement(uint i, uint j) const;
		/*!< \brief access to element in matrix at row j, column j								*/

#ifdef __NVCC__
__host__ __device__
#endif
	void loadIdentity();
		/*!< \brief load identity matrix														*/

#ifdef __NVCC__
__host__ __device__
#endif
	void loadZeroes();
		/*!< \brief initialize matrix with zeros												*/

#ifdef __NVCC__
__host__ __device__
#endif
	void loadTransMat(const Vec<rows-1, T> &vec);
		/*!< \brief creates a translation matrix in homogeneous coordinates 					*/

#ifdef __NVCC__
__host__ __device__
#endif
	void switchRows(uint i, uint j);
		/*!< \brief changes rows i and j														*/

#ifdef __NVCC__
__host__ __device__
#endif
	void addRow(uint dest, uint source, T scale);
		/*! \brief add scale * row(source) to row(dest)											*/

#ifdef __NVCC__
__host__ __device__
#endif
	void scale(T scale);
		/*!< \brief scales every entry															*/

#ifdef __NVCC__
__host__ __device__
#endif
	void scaleRow(uint row, T scale);
		/*!< \brief scales a specific row by scale												*/

#ifdef __NVCC__
__host__ __device__
#endif
	uint solveSLE(Vec<rows, hcFloat> &result);
		/*!< \brief solves a system of linear equations	using gauss' method						*/

#ifdef __NVCC__
__host__ __device__
#endif
	int translate(const Vec<rows-1, T> &vec);
		/*!< \brief assumes homogeneous TF-matrix, translates scene by vec						*/

#ifdef __NVCC__
__host__ __device__
#endif
	Matrix<cols, rows, T> transpose();
		/*!< \brief transposes matrix															*/

#ifdef __NVCC__
__host__ __device__
#endif
	void dump() const;

#ifdef __NVCC__
__host__ __device__
#endif
	void subdump();
};
//---------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------
//---			Matrix end
//---------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------

template<uint rows1, uint cols, uint cols2, class T>
#ifdef __NVCC__
__host__ __device__
#endif
Matrix<rows1, cols2, T> operator*(const Matrix<rows1, cols, T> &lhs, const Matrix<cols, cols2, T> &rhs);

template<uint rows, uint cols, class T>
#ifdef __NVCC__
__host__ __device__
#endif
Matrix<rows, cols, T> operator*(Matrix<rows, cols, T> lhs, float rhs);

template<uint rows, uint cols, class T>
#ifdef __NVCC__
__host__ __device__
#endif
Matrix<rows, cols, T> operator/(Matrix<rows, cols, T> lhs, float rhs);

template<uint rows, uint cols, class T>
#ifdef __NVCC__
__host__ __device__
#endif
Matrix<rows, cols, T> operator*(float lhs, Matrix<rows, cols, T> rhs);

template<uint rows, uint cols, class T>
#ifdef __NVCC__
__host__ __device__
#endif
Matrix<rows, cols, T> operator/(float lhs, Matrix<rows, cols, T> rhs);

template<uint rows, uint cols, class T>
#ifdef __NVCC__
__host__ __device__
#endif
Matrix<rows, cols, T> operator*(Matrix<rows, cols, T> lhs, double rhs);

template<uint rows, uint cols, class T>
#ifdef __NVCC__
__host__ __device__
#endif
Matrix<rows, cols, T> operator/(Matrix<rows, cols, T> lhs, double rhs);

template<uint rows, uint cols, class T>
#ifdef __NVCC__
__host__ __device__
#endif
Matrix<rows, cols, T> operator*(double lhs, Matrix<rows, cols, T> rhs);

template<uint rows, uint cols, class T>
#ifdef __NVCC__
__host__ __device__
#endif
Matrix<rows, cols, T> operator/(double lhs, Matrix<rows, cols, T> rhs);

template<uint rows, uint cols, class T>
#ifdef __NVCC__
__host__ __device__
#endif
Matrix<rows, cols, T> operator*(Matrix<rows, cols, T> lhs, long double rhs);

template<uint rows, uint cols, class T>
#ifdef __NVCC__
__host__ __device__
#endif
Matrix<rows, cols, T> operator/(Matrix<rows, cols, T> lhs, long double rhs);

template<uint rows, uint cols, class T>
#ifdef __NVCC__
__host__ __device__
#endif
Matrix<rows, cols, T> operator*(long double lhs, Matrix<rows, cols, T> rhs);

template<uint rows, uint cols, class T>
#ifdef __NVCC__
__host__ __device__
#endif
Matrix<rows, cols, T> operator/(long double lhs, Matrix<rows, cols, T> rhs);

template<uint rows, uint cols, class T>
#ifdef __NVCC__
__host__ __device__
#endif
Matrix<rows, cols, T> operator+(Matrix<rows, cols, T> lhs, const Matrix<rows, cols, T> &rhs);

template<uint rows, uint cols, class T>
#ifdef __NVCC__
__host__ __device__
#endif
Matrix<rows, cols, T> operator-(Matrix<rows, cols, T> lhs, const Matrix<rows, cols, T> &rhs);

template<uint rows, uint cols, class T>
Matrix<rows, cols, T>::Matrix()
{
    for(uint i=0; i<rows*cols; ++i)
        content[i] = 0.0;
}

template<uint rows, uint cols, class T>
Matrix<rows, cols, T>::Matrix(const Vec<rows, T> &v1, const Vec<rows, T> &v2)
{
	for(uint i = 0; i<rows; ++i)
	{
		content[i*cols] 	= v1[i];
		content[i*cols+1] 	= v2[i];
	}

}

template<uint rows, uint cols, class T>
Matrix<rows, cols, T>::Matrix(const Vec<rows, T> &v1, const Vec<rows, T> &v2, const Vec<rows, T> &v3)
{
	for(uint i = 0; i<rows; ++i)
	{
		content[i*cols] 	= v1[i];
		content[i*cols+1] 	= v2[i];
		content[i*cols+2] 	= v3[i];
	}
}

template<uint rows, uint cols, class T>
template<class S>
Matrix<rows, cols, T>::Matrix(const Matrix<rows, cols, S> &other)
{
	*this = other;
}

template<uint rows, uint cols, class T>
Matrix<rows, cols, T>::~Matrix(){}

template<uint rows, uint cols, class T>
template<class S>
Matrix<rows, cols, T> &Matrix<rows, cols, T>::operator=(const Matrix<rows, cols, S> &other){

	//if(this == &other)
	//	return *this;

	for(uint i = 0; i<rows*cols; ++i)
		content[i] = other.content[i];

	return *this;
}

template<uint rows, uint cols, class T>
Matrix<rows, cols, T> &Matrix<rows, cols, T>::operator*=(T scale)
{
    for(uint i=0; i<rows; ++i)
        for(uint j=0; j<cols; ++j)
            this->operator()(i, j) *= scale;
    return *this;
}

template<uint rows, uint cols, class T>
Matrix<rows, cols, T> &Matrix<rows, cols, T>::operator/=(T scale)
{
    for(uint i=0; i<rows; ++i)
        for(uint j=0; j<cols; ++j)
            this->operator()(i, j) /= scale;
    return *this;
}

template<uint rows, uint cols, class T>
Matrix<rows, cols, T> &Matrix<rows, cols, T>::operator+=(const Matrix<rows, cols, T> &other)
{
    Matrix &lhs = *this;

    for(uint i=0; i<rows; ++i)
        for(uint j=0; j<cols; ++j)
            lhs(i, j) += other[i * cols + j];

    return *this;
}

template<uint rows, uint cols, class T>
Matrix<rows, cols, T> &Matrix<rows, cols, T>::operator-=(const Matrix<rows, cols, T> &other)
{
    Matrix &lhs = *this;

    for(uint i=0; i<rows; ++i)
        for(uint j=0; j<cols; ++j)
            lhs(i, j) -= other[i * cols + j];
    return *this;
}

template<uint rows, uint cols, class T>
Vec<rows, T> Matrix<rows, cols, T>::operator*(const Vec<cols, T> &vec)
{
    Vec<rows, T> retval;
    for(uint i=0; i<rows; ++i)
        for(uint j=0; j<cols; ++j)
            retval(i) += this->operator()(i, j) * vec.content[j];
    return retval;
}

template<uint rows, uint cols, class T>
T &Matrix<rows, cols, T>::operator()(uint i, uint j)
{
	return(content[i*cols + j]);
}

template<uint rows, uint cols, class T>
T Matrix<rows, cols, T>::operator[](uint n) const
{
	return(content[n]);
}

template<uint rows, uint cols, class T>
T Matrix<rows, cols, T>::getElement(uint i, uint j) const
{
	return(content[i*cols + j]);
}

template<uint rows, uint cols, class T>
void Matrix<rows, cols, T>::loadIdentity()
{
    if (rows!=cols)
        printf("ERROR! MATRIX::loadIdentity: Row and column number not equal!\n");
    else
    {
        for(uint i=0;i<rows;++i)
            for(uint j=0;j<cols;++j)
            	content[i*cols+j] = i==j ? 1.0 : 0.0;
    }
}

template<uint rows, uint cols, class T>
void Matrix<rows, cols, T>::loadZeroes()
{
    for(uint i=0;i<rows;++i)
        for(uint j=0;j<cols;++j)
            content[i*cols+j] = 0.0;
}

template<uint rows, uint cols, class T>
void Matrix<rows, cols, T>::loadTransMat(const Vec<rows-1, T> &vec)
{
	this->loadIdentity();
	for(uint i=0;i<rows-1;++i)
		content[(i+1)*cols-1] = vec.content[i];
}

template<uint rows, uint cols, class T>
void Matrix<rows, cols, T>::switchRows(uint i, uint j)
{
	for(uint k=0; k<cols; ++k)
	{
		T temp 			= content[i*cols + k];
		content[i*cols + k] = content[j*cols + k];
		content[j*cols + k] = temp;
	}
}

template<uint rows, uint cols, class T>
void Matrix<rows, cols, T>::addRow(uint dest, uint source, T scale)
{
	for (uint i=0; i<cols; ++i)
		content[dest*cols + i] += scale * content[source*cols+i];
}


/*!	Solves a System of Linear Equations (Ax=y) posed by the matrix "Ay" for x, where the rightmost column is
 * 	interpreted as the solution vector y. The Gaussian elimination algorithm is used.
 * 	So far, only unique solutions are stored in the vector result. The matrix Ay is not altered but copied to a local
 * 	object
 *
 * 	\retval 0 there is no solution or ill posed problem
 * 	\retval 1 unique solution, stored now in result
 * 	\retval 2 there is an infinte number of solution
 *
 * 	only in case 1 does the result vector contain useful data
 */
template<uint rows, uint cols, class T>
uint Matrix<rows, cols, T>::solveSLE(Vec<rows, hcFloat> &result)
{
	//hcFloat numeps = 1E-20;
	hcFloat numeps = 1E-6;

	Matrix A(*this);
	//printf("Starting SLE:\n");
	//A.dump();

	uint i = 0;
	uint j = 0;
	while (i<rows)                    		// create triangular matrix
	{
		if (fabs(A(i, i)) < numeps)		// pivot element to small, has to be substituted
		{
			//printf("SolveSLE: element too small...\n");
			j = i+1;
			while((j < rows) && (fabs(A(j, i)) < numeps))
				++j;

			if (j < rows)                   // find some row that can be exchanged with current row
				A.switchRows(i,j);
			else                           	// there is no such row, check if SLE is unsolvable
			{
				j = i;
				while(j<cols-1 && fabs(A(i, j)) < numeps)
					++j;

				if(j < cols-1)              // ok, there is no row to be exchanged, but maybe SLE is still solvable
					++i;
				else
				{
					//if (fabs(A[(i+1)*cols - 1]) < numeps)
					if (fabs(A(i, cols-1) < numeps))
						++i;
					else                    // SLE not solvable
						return 0;
				}
			}
		}
		else
		{
			j=i+1;							// Gauss step adding a scaled row to all the others that have non-0 elements in column i
			while(j<rows)
			{
				if(fabs(A(j, i)) > numeps)
					A.addRow(j, i, -A(j,i) / A(i,i));
				++j;
			}
			++i;
		}
	}
	//printf("triang:\n");
	//A.dump();

	i=rows-1;
	while(true)
	{
		if(fabs(A(i, i)) < numeps)
			return 2;						// SLE has infinitely many solutions
		else
		{
			//result[i] = A[(i+1)*cols - 1];
			result(i) = A(i, cols-1);

			for(j=rows-1; j>i; --j)
			{
				result(i) -= result[j] * A(i, j);
				//printf("%u %E\n", i, result[i]);
			}

			//result[i] /= A[i*(cols+1)];
			result(i) /= A(i, i);

			//printf("final: %E\n", result[i]);

			if (i == 0)
				return 1;					// unique solution found and stored in result

			--i;
		}
	}
}

template<uint rows, uint cols, class T>
void Matrix<rows, cols, T>::scale(T scale)
{
    for(uint i=0;i<rows;++i)
        for(uint j=0;j<cols;++j)
            content[i*cols+j] *= scale;
}

template<uint rows, uint cols, class T>
void Matrix<rows, cols, T>::scaleRow(uint row, T scale)
{
	if(row >= rows)
	{
		printf("ERROR! hcMatrix::scaleRow: Requested row (%u) outside scope (#rows: %u)\n", row, rows);
		return;
	}

    for(uint j=0;j<cols;++j)
        content[row * cols + j] *= scale;
}

template<uint rows, uint cols, class T>
int Matrix<rows, cols, T>::translate(const Vec<rows-1, T> &vec)
{
    Matrix<rows, cols, T> temp;
    temp.loadIdentity();

    for(uint i=0; i<rows-1; ++i)
        temp.content[(i+1)*cols - 1] = vec.content[i];

    *this = temp * *this;
    return 1;
}

template<uint rows, uint cols, class T>
Matrix<cols, rows, T> Matrix<rows, cols, T>::transpose()
{
    Matrix<cols, rows, T> retval;

    for(uint i=0; i<rows; ++i)
    	for(uint j=0; j<cols; ++j)
    		retval(j,i) = this->operator ()(i, j);
    return retval;
}

template<uint rows, uint cols, class T>
void Matrix<rows, cols, T>::dump() const
{
    printf("Dumping (%u x %u)-hcMatrix:\n", rows, cols);
	for(uint i=0; i<rows; ++i)
	{
		for(uint j=0; j<cols; ++j)
			printf("%E\t",content[i*cols+j]);
		printf("\n");
	}
}

template<uint rows, uint cols, class T>
void Matrix<rows, cols, T>::subdump()
{
    uint numEl 		= 4;
    uint lastCol 	= 0;
    uint i;

    while(lastCol < cols)
    {
        printf("Here comes a submatrix:\n");
        i=0;
        while(i<rows)
        {
            for(uint j=lastCol; j<lastCol+numEl; ++j)
                printf("%E\t", this->content[i*cols+j]);
            printf("\n");
            ++i;
        }
        printf("\n");
        lastCol = lastCol + numEl;
    }
}

template<uint rows, uint cols, class T>
Matrix<rows, cols, T> operator*(float lhs, Matrix<rows, cols, T> rhs)
{
    rhs *= lhs;
    return rhs;
}

template<uint rows, uint cols, class T>
Matrix<rows, cols, T> operator*(Matrix<rows, cols, T> lhs, float rhs)
{
    lhs *= rhs;
    return lhs;
}

template<uint rows, uint cols, class T>
Matrix<rows, cols, T> operator*(double lhs, Matrix<rows, cols, T> rhs)
{
    rhs *= lhs;
    return rhs;
}

template<uint rows, uint cols, class T>
Matrix<rows, cols, T> operator*(Matrix<rows, cols, T> lhs, double rhs)
{
    lhs *= rhs;
    return lhs;
}

template<uint rows, uint cols, class T>
Matrix<rows, cols, T> operator*(long double lhs, Matrix<rows, cols, T> rhs)
{
    rhs *= lhs;
    return rhs;
}

template<uint rows, uint cols, class T>
Matrix<rows, cols, T> operator*(Matrix<rows, cols, T> lhs, long double rhs)
{
    lhs *= rhs;
    return lhs;
}

template<uint rows, uint cols, class T>
Matrix<rows, cols, T> operator/(float lhs, Matrix<rows, cols, T> rhs)
{
    rhs /= lhs;
    return rhs;
}

template<uint rows, uint cols, class T>
Matrix<rows, cols, T> operator/(Matrix<rows, cols, T> lhs, float rhs)
{
    lhs /= rhs;
    return lhs;
}

template<uint rows, uint cols, class T>
Matrix<rows, cols, T> operator/(double lhs, Matrix<rows, cols, T> rhs)
{
    rhs /= lhs;
    return rhs;
}

template<uint rows, uint cols, class T>
Matrix<rows, cols, T> operator/(Matrix<rows, cols, T> lhs, double rhs)
{
    lhs /= rhs;
    return lhs;
}

template<uint rows, uint cols, class T>
Matrix<rows, cols, T> operator+(Matrix<rows, cols, T> lhs, const Matrix<rows, cols, T> &rhs){

    lhs += rhs;
    return lhs;
}

template<uint rows, uint cols, class T>
Matrix<rows, cols, T> operator-(Matrix<rows, cols, T> lhs, const Matrix<rows, cols, T> &rhs){

    lhs -= rhs;
    return lhs;
}

template<uint rows1, uint cols, uint cols2, class T>
Matrix<rows1, cols2, T> operator*(const Matrix<rows1, cols, T> &lhs, const Matrix<cols, cols2, T> &rhs){

	Matrix<rows1, cols2, T> retval;

	for(uint i=0; i<rows1; ++i)
		for(uint j=0; j<cols2; ++j)
			for(uint k=0; k<cols; ++k)
				retval(i, j) += lhs.getElement(i, k) * rhs.getElement(k, j);

	return retval;
}

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          MatrixNxN
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------


/*! \brief implementation of a square matrix
 *
 */
template<uint dim, class T>
class MatrixNxN : public Matrix<dim, dim, T>{
public:

#ifdef __NVCC__
__host__ __device__
#endif
	MatrixNxN();
		/*!< \brief std constructor																*/

template<class S>
#ifdef __NVCC__
__host__ __device__
#endif
	MatrixNxN(const MatrixNxN<dim, S> &other);
		/*!< \brief cpy constructor																*/

#ifdef __NVCC__
__host__ __device__
#endif
	~MatrixNxN();
		/*!< \brief destructor																	*/

#ifdef __NVCC__
__host__ __device__
#endif
	MatrixNxN<dim, T> &operator=(const MatrixNxN<dim, T> &other);
		/*!< \brief assignment operator															*/

#ifdef __NVCC__
__host__ __device__
#endif
	double det();
		/*!< \brief computes determinant of matrix (double precision)							*/

#ifdef __NVCC__
__host__ __device__
#endif
	int invert();
		/*!< \brief computes inverse of matrix													*/

#ifdef __NVCC__
__host__ __device__
#endif
	void loadIdentity();
		/*!< \brief initializes matrix with ones on main diagonal								*/

#ifdef __NVCC__
__host__ __device__
#endif
    void scaleAllDim(T scale);
        /*!< \brief applies a scale operation on all dimensions                                 */

#ifdef __NVCC__
__host__ __device__
#endif
	void scaleHom(T scale);
		/*!< \brief assumes homogeneous TF-matrix, scales the scene isotropically				*/

#ifdef __NVCC__
__host__ __device__
#endif
	void scaleAxis(uint n, T scale);
		/*!< \brief assumes homogeneous TF-matrix, scales along a specific axis					*/

#ifdef __NVCC__
__host__ __device__
#endif
	void translate(const Vec<dim-1, T> &v);
		/*!< \brief assumes homogeneous TF-matrix, translates along vec v in N-1 dimensions		*/

#ifdef __NVCC__
__host__ __device__
#endif
    void dump() const;
};

template<uint dim, class T, class S>
#ifdef __NVCC__
__host__ __device__
#endif
MatrixNxN<dim, T> operator*(const MatrixNxN<dim, T> &lhs, const MatrixNxN<dim, S> &rhs);


template<uint dim, class T>
MatrixNxN<dim, T>::MatrixNxN() : Matrix<dim,dim, T>(){}

template<uint dim, class T>
template<class S>
MatrixNxN<dim, T>::MatrixNxN(const MatrixNxN<dim, S> &other) : Matrix<dim, dim, T>(other){}

template<uint dim, class T>
MatrixNxN<dim, T>::~MatrixNxN(){}

template<uint dim, class T>
MatrixNxN<dim, T> &MatrixNxN<dim, T>::operator=(const MatrixNxN<dim, T> &other)
{
    if(this == &other)	return *this;

	Matrix<dim, dim, T>::operator=(other);

	return *this;
}

template<uint dim, class T>
double MatrixNxN<dim, T>::det()
{
    uint j 			= 0; // TODO: what does this do?
    int sign		= 0;
    double result 	= 0.0;

    if(dim == 2)	return((double)this->content[0] * (double)this->content[3] - (double)this->content[2] * (double)this->content[1]);
    else
    {
		for(uint i=0;i<dim;++i)
		{
			if((i + 1 + j +1) % 2 == 0)
				sign = 1;
			else
				sign = -1;

			uint c = 0;
			MatrixNxN<dim-(dim>2 ? 1 : 0), T> tempmat;
			for(uint k=0;k<dim;++k)
				for(uint l=0;l<dim;++l)
					if((k!=i) && (l!=j))
						tempmat.content[c++] = this->content[k*dim+l];

			double det = tempmat.det();
			result += sign * (double)this->content[i*dim] * det;
		}
		return result;
    }
}


template<uint dim, class T>
int MatrixNxN<dim, T>::invert()
{
	Matrix<dim, 2*dim, T> workhorse;

	for(uint i=0;i<dim;++i)
		for(uint j=0;j<dim;++j)
			workhorse.content[i*(2*dim) + j] = this->content[i*dim+j];

	for(uint i=0;i<dim;++i)
		workhorse.content[i*(2*dim) + dim + i] = 1.0;

	uint j = 0;
	uint i = 0;
	while (i<dim)                            // create triangular matrix
	{
		if (fabs(workhorse.content[i*(2*dim) + i]) < num_eps)         // there is some diagonal element which is zero
		{
			j = i+1;
			while((j < dim) && (fabs(workhorse.content[j*(2*dim) + i]) < num_eps))
				++j;

			if (j < dim)                     // find some row that can be exchanged with current row
				workhorse.switchRows(i,j);
			else
			{
				printf("WARNING! MatrixNxN::invert: No inverse found!\n");
				return 0;
			}
		}
		else
		{
			j=0;
			while(j<dim)
			{
				if(j!=i)
				{
					if(fabs(workhorse.content[j*(2*dim) + i]) > num_eps)
						workhorse.addRow(j, i, -workhorse.content[j*(2*dim)+i] / workhorse.content[i*(2*dim)+i]);
				}
				else
					workhorse.scaleRow(i, 1 / workhorse.content[i*(2*dim)+i]);
				++j;
			}
			++i;
		}
	}

	for(i=0;i<dim;++i)
		for(j=0;j<dim;++j)
			this->content[i*dim+j] = workhorse.content[i*(2*dim) + dim + j];
	return 1;
}

template<uint dim, class T>
void MatrixNxN<dim, T>::loadIdentity(){

    for(uint i=0;i<dim;++i)
        for(uint j=0;j<dim;++j)
        {
            if(i==j)

                this->content[i*dim+j] = 1.0;
            else
                this->content[i*dim+j] = 0.0;
        }

}

template<uint dim, class T>
void MatrixNxN<dim, T>::scaleAxis(uint n, T scale){

	if(n > dim)
        printf("ERROR! MatrixNxN::scaleAxis: Axis to be scaled (%u) does not exist (dim=%u)!\n", n, dim);
    else
        for(uint j=0; j<dim; ++j)
            this->content[n*dim+j] *= scale;
}

template<uint dim, class T>
void MatrixNxN<dim, T>::scaleAllDim(T scale){

    MatrixNxN<dim, T> scaleMat;
    scaleMat.loadIdentity();
    for(uint i=0; i<dim-1; ++i)
        scaleMat(i, i) *= scale;

    *this = scaleMat * *this;
}

template<uint dim, class T>
void MatrixNxN<dim, T>::scaleHom(T scale){
    this->content[dim*dim-1] /= scale;
}

template<uint dim, class T>
void MatrixNxN<dim, T>::translate(const Vec<dim-1, T> &v){

	for(uint i=1; i<=dim-1; ++i)
		this->content[i*dim-1] += v.content[i-1];
}

template<uint dim, class T>
void MatrixNxN<dim, T>::dump() const{

    printf("Dumping MatrixNxN with dim %u:\n", dim);
    Matrix<dim, dim, T>::dump();
}

template<uint dim, class T, class S>
MatrixNxN<dim, T> operator*(const MatrixNxN<dim, T> &lhs, const MatrixNxN<dim, S> &rhs){

	MatrixNxN<dim, T> retval;

	for(uint i=0; i<dim; ++i)
		for(uint j=0; j<dim; ++j)
			for(uint k=0; k<dim; ++k)
				retval(i, j) += lhs.getElement(i, k) * rhs.getElement(k, j);

	return retval;
}


class Matrix2x2 : public MatrixNxN<2, hcFloat>{
public:

#ifdef __NVCC__
__host__ __device__
#endif
    Matrix2x2();							/*!< \brief std constructor				*/

#ifdef __NVCC__
__host__ __device__
#endif
    Matrix2x2(const Matrix2x2 &other);		/*!< \brief cpy constructor				*/

template<class S>
#ifdef __NVCC__
__host__ __device__
#endif
    Matrix2x2(const MatrixNxN<2, S> &other);		/*!< \brief pseudo cpy constructor				*/

#ifdef __NVCC__
__host__ __device__
#endif
	~Matrix2x2();

#ifdef __NVCC__
__host__ __device__
#endif
    Matrix2x2 &operator=(const Matrix2x2 &other);

template<class S>
#ifdef __NVCC__
__host__ __device__
#endif
    Matrix2x2 &operator=(const MatrixNxN<2, S> &other);

template<class S>
#ifdef __NVCC__
__host__ __device__
#endif
    Matrix2x2 &operator=(const Matrix<2, 2, S> &other);

template<class S>
#ifdef __NVCC__
__host__ __device__
#endif
    void scalex(S scale);

template<class S>
#ifdef __NVCC__
__host__ __device__
#endif
    void scaley(S scale);

template<class S>
#ifdef __NVCC__
__host__ __device__
#endif
    void scale(S scalex, S scaley);

};

/*! \brief implements a 3x3 matrix for use, e.g., with homogeneus 2D space
 *
 */
class Matrix3x3 : public MatrixNxN<3, hcFloat>{
public:

#ifdef __NVCC__
__host__ __device__
#endif
    Matrix3x3();

#ifdef __NVCC__
__host__ __device__
#endif
    Matrix3x3(const Matrix3x3 &other);

template<class S>
#ifdef __NVCC__
__host__ __device__
#endif
    Matrix3x3(const MatrixNxN<3, S> &other);

template<class S>
#ifdef __NVCC__
__host__ __device__
#endif
    Matrix3x3(const Vec<3, S> &vec1, const Vec<3, S> &vec2, const Vec<3, S> &vec3);

#ifdef __NVCC__
__host__ __device__
#endif
    ~Matrix3x3();

#ifdef __NVCC__
__host__ __device__
#endif
    Matrix3x3 &operator=(const Matrix3x3 &other);

template<class S>
#ifdef __NVCC__
__host__ __device__
#endif
    Matrix3x3 &operator=(const MatrixNxN<3, S> &other);

#ifdef __NVCC__
__host__ __device__
#endif
    void convertSphericalToCartesian(const Vec3D &cartPos);

#ifdef __NVCC__
__host__ __device__
#endif
    void convertCartesianToSpherical(const Vec3D &cartPos);
    	/*! \brief loads the conversion matrix to convert from cartesian coords to spherical coords */

template<class S>
#ifdef __NVCC__
__host__ __device__
#endif
    void scalex(S scale);

template<class S>
#ifdef __NVCC__
__host__ __device__
#endif
    void scaley(S scale);

template<class S>
#ifdef __NVCC__
__host__ __device__
#endif
    void scalez(S scale);

template<class S>
#ifdef __NVCC__
__host__ __device__
#endif
    void scale(S scalex, S scaley, S scalez);

template<class S>
#ifdef __NVCC__
__host__ __device__
#endif
    void loadRotationX(S theta);

template<class S>
#ifdef __NVCC__
__host__ __device__
#endif
    void loadRotationY(S theta);

template<class S>
#ifdef __NVCC__
__host__ __device__
#endif
    void loadRotationZ(S theta);

template<class S>
#ifdef __NVCC__
__host__ __device__
#endif
	void rotateAroundAxis(Vec3D axis, S angle);

template<class S>
#ifdef __NVCC__
__host__ __device__
#endif
    void loadEulerTransform(S omega, S theta, S phi);	// = Rz(phi)*Rx(theta)*Rz(omega) (see Heliospheric coordinate systems, Fr��nz, Harper 2002)

};


class hcQuaternion;
/*! \brief implements a 4x4 matrix for use, e.g., with homogeneous 3D space
 *
 */
class Matrix4x4 : public MatrixNxN<4, hcFloat>{
public:

#ifdef __NVCC__
__host__ __device__
#endif
    Matrix4x4();									/*!< \brief std constructor					*/

#ifdef __NVCC__
__host__ __device__
#endif
    Matrix4x4(const Matrix4x4 &other);				/*!< \brief cpy constructor					*/

template<class S>
#ifdef __NVCC__
__host__ __device__
#endif
    Matrix4x4(const MatrixNxN<4, S> &other);			/*!< \brief pseudo-cpy constructor			*/

#ifdef __NVCC__
__host__ __device__
#endif
    ~Matrix4x4();									/*!< \brief std constructor					*/

#ifdef __NVCC__
__host__ __device__
#endif
	Matrix4x4 &operator=(const Matrix4x4 &other);

template<class S>
#ifdef __NVCC__
__host__ __device__
#endif
	Matrix4x4 &operator=(const MatrixNxN<4, S> &other);

using MatrixNxN::scale;
using MatrixNxN::translate;

#ifdef GUI
#ifdef __NVCC__
__host__ __device__
#endif
    void set(const aiMatrix4x4 &mat);
#endif

#ifdef __NVCC__
__host__ __device__
#endif
    void loadTFMatrix(const Vec3D &right, const Vec3D &up, const Vec3D &back, const Vec3D &pos);

template<class S>
#ifdef __NVCC__
__host__ __device__
#endif
    void scalex(S scale);

template<class S>
#ifdef __NVCC__
__host__ __device__
#endif
    void scaley(S scale);

template<class S>
#ifdef __NVCC__
__host__ __device__
#endif
    void scalez(S scale);

template<class S>
#ifdef __NVCC__
__host__ __device__
#endif
    void scalew(S scale);

template<class S>
#ifdef __NVCC__
__host__ __device__
#endif
    void scale(S scalex, S scaley, S scalez, S scalew);

template<class S>
#ifdef __NVCC__
__host__ __device__
#endif
    void rotatex(S phi);

template<class S>
#ifdef __NVCC__
__host__ __device__
#endif
    void rotatey(S phi);

template<class S>
#ifdef __NVCC__
__host__ __device__
#endif
    void rotatez(S phi);

#ifdef __NVCC__
__host__ __device__
#endif
    void rotateAtoB(const Vec3D a, const Vec3D b);
        /*!< \brief creates a rotation matrix that rotates a onto b                 */

template<class S>
#ifdef __NVCC__
__host__ __device__
#endif
    void translate(S x, S y, S z);

#ifdef GUI
#ifdef __NVCC__
__host__ __device__
#endif
    void scaleRotTrans(const Vec3D &scale, const hcQuaternion &rot, const Vec3D &trans);
#endif
};


/*! \brief implements a 5x5 matrix for use, e.g., with homogeneous 4D space
 *
 * 	TODO: not tested
 *
 */
class Matrix5x5 : public MatrixNxN<5, hcFloat>{
public:

#ifdef __NVCC__
__host__ __device__
#endif
    Matrix5x5();

#ifdef __NVCC__
__host__ __device__
#endif
    Matrix5x5(const Matrix5x5 &other);

template<class S>
#ifdef __NVCC__
__host__ __device__
#endif
    Matrix5x5(const MatrixNxN<5, S> &other);

#ifdef __NVCC__
__host__ __device__
#endif
    ~Matrix5x5();

#ifdef __NVCC__
__host__ __device__
#endif
	Matrix5x5 &operator=(const Matrix5x5 &other);

template<class S>
#ifdef __NVCC__
__host__ __device__
#endif
	Matrix5x5 &operator=(const MatrixNxN<5, S> &other);

#ifdef __NVCC__
__host__ __device__
#endif
    void loadTFMatrix(const Vec4D &right, const Vec4D &up, const Vec4D &back, const Vec4D &over, const Vec4D &pos);

template<class S>
#ifdef __NVCC__
__host__ __device__
#endif
    void rotateXY(S phi);

template<class S>
#ifdef __NVCC__
__host__ __device__
#endif
    void rotateXZ(S phi);

template<class S>
#ifdef __NVCC__
__host__ __device__
#endif
    void rotateYZ(S phi);

template<class S>
#ifdef __NVCC__
__host__ __device__
#endif
    void rotateXW(S phi);

template<class S>
#ifdef __NVCC__
__host__ __device__
#endif
    void rotateYW(S phi);

template<class S>
#ifdef __NVCC__
__host__ __device__
#endif
    void rotateZW(S phi);

template<class S>
#ifdef __NVCC__
__host__ __device__
#endif
    void translate(S x, S y, S z);
};

#endif

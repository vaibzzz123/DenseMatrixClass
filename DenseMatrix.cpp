#include <iostream>
#include <vector>
#include <cmath>


static double DetHelper(std::vector <std::vector<double>> A, int Rows, int Cols);
static double Cofactor(std::vector<std::vector<double>> A, int Row, int Col);

//TODO: have consistent numbering scheme (1->n or 0->n-1)
//TODO: make public variables private and make getter methods for them;

class DenseMatrix {
public:
	int Rows;
	int Cols;
	std::vector<std::vector<double>> Matrix;

	DenseMatrix(int R, int C);

	//construct matrix from vector of vectors
	//requires all vectors to be of equal size
	static DenseMatrix *ConstructRow(std::vector<std::vector<double>>* A);
	static DenseMatrix *ConstructCol(std::vector<std::vector<double>>* A);

	//return a row/column as a C++ vector
	bool IsValidRow(DenseMatrix *M, int Row);
	static std::vector<double>* ReturnRow(DenseMatrix *M, int Row);
	bool IsValidCol(DenseMatrix *M, int Col);
	static std::vector<double>* ReturnCol(DenseMatrix *M, int Col);

	//Elementary Row Operations
	static void RowIPlusCRowJ(DenseMatrix* M, int i, int j, int c);
	static void RowISwapRowJ(DenseMatrix *M, int i, int j);
	static void CRowI(DenseMatrix *M, int i, int c);

	//construct matrix columns/rows from vectors

	//linear combinations of matrices
	//returns pointer to new dynamic memory, caller must delete
	static DenseMatrix* ScalarMultiplication(DenseMatrix *M, double c); //why  is static required? It just works when I do it???
	bool IsAddable(DenseMatrix *A, DenseMatrix *B);
	//returns pointer to new dynamic memory, caller must delete
	static DenseMatrix* MatrixAddition(DenseMatrix *A, DenseMatrix *B);
	//makes n x n identity matrix, returns pointer to new dynamic memory, caller must delete
	static DenseMatrix* IdentityMatrix(int n);

	//matrix vector multiplication
	bool IsMVMultipliable(DenseMatrix *M, std::vector<double> *V);
	static std::vector<double>* MatrixVectorMultiplication(DenseMatrix *M, std::vector<double> *V);

	//matrix matrix multiplication
	bool IsMMMultipliable(DenseMatrix *A, DenseMatrix *B);
	static DenseMatrix* MatrixMatrixMultiplication(DenseMatrix *A, DenseMatrix *B);

	//determinant of a matrix
	static double Determinant (DenseMatrix *M);

	//inverse of a matrix
	bool IsInvertible(DenseMatrix *M);
	static DenseMatrix* MatrixInverse(DenseMatrix *M);

	//some IO functions
	static void PrintMatrix(DenseMatrix *M);
	static void PrintVector(std::vector<double>* V);

private:

};

DenseMatrix::DenseMatrix(int R, int C)
{
	Matrix = std::vector<std::vector<double>>(R, std::vector<double>(C, 0)); //initializes all values to 0 (HOW DOES THIS EVEN WORK???)

	this->Rows = R;
	this->Cols = C;
}

DenseMatrix * DenseMatrix::ConstructRow(std::vector<std::vector<double>>* A)
{
	int Rows = A->size();
	int Cols = (*A)[0].size();
	DenseMatrix *M = new DenseMatrix(Rows, Cols);
	for (int i = 0; i < Rows; ++i)
	{
		for (int j = 0; j < Cols; ++j)
		{
			M->Matrix[i][j] = (*A)[i][j];
		}
	}
	return M;
}

DenseMatrix * DenseMatrix::ConstructCol(std::vector<std::vector<double>>* A)
{
	int Rows = (*A)[0].size();
	int Cols = A->size();
	DenseMatrix *M = new DenseMatrix(Rows, Cols);
	for (int i = 0; i < Rows; ++i)
	{
		for (int j = 0; j < Cols; ++j)
		{
			M->Matrix[i][j] = (*A)[j][i];
		}
	}
	return M;
}

bool DenseMatrix::IsValidRow(DenseMatrix * M, int Row)
{
	return (0 <= Row) && (Row < M->Rows);
}

std::vector<double>* DenseMatrix::ReturnRow(DenseMatrix * M, int Row)
{
	std::vector<double> *N = new std::vector<double>(M->Cols, 0.0);
	int Cols = M->Cols;
	for (int i = 0; i < Cols; ++i)
	{
		(*N)[i] = M->Matrix[Row][i];
	}
	return N;
}

bool DenseMatrix::IsValidCol(DenseMatrix * M, int Col)
{
	return (0 <= Col) && (Col < M->Cols);
}

std::vector<double>* DenseMatrix::ReturnCol(DenseMatrix * M, int Col)
{
	std::vector<double> *N = new std::vector<double>(M->Rows, 0.0);
	int Rows = M->Rows;
	for (int i = 0; i < Rows; ++i)
	{
		(*N)[i] = M->Matrix[i][Col];
	}
	return N;
}

void DenseMatrix::RowIPlusCRowJ(DenseMatrix * M, int i, int j, int c)
{
	int Cols = M->Cols;
	for (int t = 0; t < Cols; ++t)
	{
		M->Matrix[i][t] = M->Matrix[i][t] + c*M->Matrix[j][t];
	}
}

void DenseMatrix::RowISwapRowJ(DenseMatrix * M, int i, int j)
{
	int Cols = M->Cols;
	std::vector<double> RowI(Cols, 0);
	for (int t = 0; t < Cols; ++t)
	{
		RowI[t] = M->Matrix[i][t]; //copying the i'th row to a vector
	}
	for (int y = 0; y < Cols; ++y)
	{
		M->Matrix[i][y] = M->Matrix[j][y]; //overwriting the i'th row with j'th row
	}
	for (int x = 0; x < Cols; ++x)
	{
		M->Matrix[j][x] = RowI[x]; //overwriting the j'th row with the i'th row copy
	}
}

void DenseMatrix::CRowI(DenseMatrix * M, int i, int c)
{
	int Cols = M->Cols;
	for (int j = 0; j < Cols; ++j)
	{
		M->Matrix[i][j] = c*M->Matrix[i][j];
	}
}

DenseMatrix* DenseMatrix::ScalarMultiplication(DenseMatrix *M, double c)
{
	int Rows = M->Rows;
	int Cols = M->Cols;
	DenseMatrix *A = new DenseMatrix(Rows, Cols);
	for (int i = 0; i < Rows; ++i)
	{
		for (int j = 0; j < Cols; ++j) {
			A->Matrix[i][j] = c*M->Matrix[i][j];
		}
	}
	return A;
}

bool DenseMatrix::IsAddable(DenseMatrix * A, DenseMatrix * B)
{
	return A->Rows == B->Rows && A->Cols == B->Cols;
}

DenseMatrix * DenseMatrix::MatrixAddition(DenseMatrix * A, DenseMatrix * B)
{
	int Rows = A->Rows;
	int Cols = A->Cols;
	DenseMatrix *C = new DenseMatrix(Rows, Cols);
	for (int i = 0; i < Rows; ++i) {
		for (int j = 0; j < Cols; ++j) {
			C->Matrix[i][j] = A->Matrix[i][j] + B->Matrix[i][j];
		}
	}
	return C;
}

DenseMatrix * DenseMatrix::IdentityMatrix(int n)
{
	int Rows = n;
	int Cols = n;
	DenseMatrix *A = new DenseMatrix(Rows, Cols);
	for (int i = 0; i < Rows; ++i)
	{
		for (int j = 0; j < Cols; ++j)
		{
			if (i == j)
			{
				A->Matrix[i][j] = 1; //diagonal elements
			}
			else //not diagonal elements
			{
				A->Matrix[i][j] = 0;
			}
		}
	}
	return A;
}

bool DenseMatrix::IsMVMultipliable(DenseMatrix * M, std::vector<double>* V)
{
	return M->Cols == V->size();
}

std::vector<double>* DenseMatrix::MatrixVectorMultiplication(DenseMatrix * M, std::vector<double>* V)
{
	//initialize new empty vector, of size M->Rows;
	int Newsize = M->Rows;
	int Rows = M->Rows;
	int Cols = M->Cols;
	std::vector<double>* B = new std::vector<double>(Newsize, 0.0);

	//for every row in the matrix
	for (int i = 0; i < Rows; ++i)
	{
		//for every column in the matrix
		for (int j = 0; j < Cols; ++j)
		{
			//the new vector's i'th element is dot product of the i'th row of matrix and vector
			double DotProduct = 0;
			for (int t = 0; t < Cols; ++t) 
			{
				DotProduct += M->Matrix[i][t] * (*V)[t];
			}
			//size_t Test = B->size();
			(*B)[i] = DotProduct;
		}
	}
	return B;
}

bool DenseMatrix::IsMMMultipliable(DenseMatrix * A, DenseMatrix * B)
{
	return A->Cols == B->Rows;
}

DenseMatrix * DenseMatrix::MatrixMatrixMultiplication(DenseMatrix * A, DenseMatrix * B)
{ //rows of A * cols of B
	DenseMatrix *C = new DenseMatrix(A->Rows, B->Cols);
	int Rows = C->Rows;
	int Cols = C->Cols;
	//for every row in the matrix
	for (int i = 0; i < Rows; ++i)
	{
		//for every column in the matrix
		for (int j = 0; j < Cols; ++j)
		{
			double DotProduct = 0;
			for (int t = 0; t < Cols; ++t)
			{
				DotProduct += A->Matrix[i][t] * B->Matrix[t][j];
			}
			C->Matrix[i][j] = DotProduct;
		}
	}
	return C;
}

double DenseMatrix::Determinant(DenseMatrix * M)
{
	return DetHelper(M->Matrix, M->Rows, M->Cols);
}

bool DenseMatrix::IsInvertible(DenseMatrix * M)
{
	return (M->Rows == M->Cols) && (DenseMatrix::Determinant(M) != 0);
}

DenseMatrix * DenseMatrix::MatrixInverse(DenseMatrix * M)
{
	int Rows = M->Rows;
	int Cols = M->Cols;
	DenseMatrix *A = new DenseMatrix(Rows, Cols); //adjugate matrix
	for (int i = 0; i < Rows; ++i)
	{
		for (int j = 0; j < Cols; ++j)
		{
			A->Matrix[i][j] = Cofactor(M->Matrix, j, i);
		}
	}
	double det = Determinant(M);
	DenseMatrix *N = ScalarMultiplication(A, 1/det);
	delete A;
	return N;
}

void DenseMatrix::PrintMatrix(DenseMatrix *M)
{
	for (auto R : M->Matrix)
	{
		for (auto x : R)
		{
			std::cout << " " << x << " ";
		}
		std::endl(std::cout);
	}
	std::endl(std::cout);
}

void DenseMatrix::PrintVector(std::vector<double>* V)
{
	int Length = static_cast<int>(V->size());
	for (int i = 0; i < Length; ++i)
	{
		std::cout << (*V)[i];
		std::endl(std::cout);
	}
	std::endl(std::cout);
}

//int main()
{
	/* DenseMatrix* M = new DenseMatrix(3,3);
	PrintMatrix(M);
	static DenseMatrix *A = DenseMatrix::ScalarMultiplication(M, 2);
	static DenseMatrix *B = DenseMatrix::ScalarMultiplication(M, 3);
	delete M;
	PrintMatrix(A);
	PrintMatrix(B);
	static DenseMatrix *C = DenseMatrix::MatrixAddition(A, B);
	PrintMatrix(C);
	std::vector<double> Odd{ 1, 3, 5 };
	PrintVector(&Odd);
	static std::vector<double>* V = DenseMatrix::MatrixVectorMultiplication(C, &Odd);
	PrintVector(V);
	static DenseMatrix* D = DenseMatrix::MatrixMatrixMultiplication(A, B);
	PrintMatrix(D);
	//std::cout << DenseMatrix::Determinant(D) << std::endl;
	delete A;
	delete B;
	delete C;
	delete V;
	static DenseMatrix *I = DenseMatrix::IdentityMatrix(3);
	DenseMatrix::PrintMatrix(I);
	std::cout << DenseMatrix::Determinant(I) << std::endl;
	DenseMatrix::RowISwapRowJ(I, 1, 2);
	//PrintMatrix(I);
	DenseMatrix::CRowI(I, 0, 2);
	//PrintMatrix(I);
	DenseMatrix::RowIPlusCRowJ(I, 1, 2, 1);
	//PrintMatrix(I);
	delete I;
	std::vector<double> A = { 5,6,7 };
	std::vector<double> B = { 3,9,8 };
	std::vector<double> C = { 4,3,2 };
	std::vector<std::vector<double>> D = { A,B,C };
	DenseMatrix *E = DenseMatrix::ConstructRow(&D);
	DenseMatrix::PrintMatrix(E);
	std::cout << DenseMatrix::Determinant(E) << std::endl;
	DenseMatrix *F = DenseMatrix::MatrixInverse(E);
	DenseMatrix::PrintMatrix(F);
	//DenseMatrix *G = DenseMatrix::MatrixMatrixMultiplication(E, F);
	//PrintMatrix(G);
	//delete G;
	delete E;
	delete F; */
	return 0;
}

double DetHelper(std::vector<std::vector<double>> A, int Rows, int Cols)
{
	double det = 0;
	if (Rows == 1) 
	{
		det = A[0][0];
	}
	else if (Rows == 2) 
	{
		det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
	}
	else 
	{
		for (int i = 0; i < Cols; ++i)
		{
			if (A[0][i] != 0)
			{
				det += A[0][i] * Cofactor(A, 0, i);
			}
		}
	}
	return det;
}

double Cofactor(std::vector<std::vector<double>> A, int Row, int Col)
{
	int R = A.size();
	int C = A[0].size();
	//let B be this new (n-1) x (n-1) 2d vector
	//initialization
	std::vector<std::vector<double>> B = std::vector<std::vector<double>>(R, std::vector<double>(C, 1));
	for (int p = 0; p < R; ++p) //copying all the elements of A into B
	{
		for (int q = 0; q < C; ++q)
		{
			B[p][q] = A[p][q];
		}
	}
	//delete row Row from B
	B.erase(B.begin() + Row);
	//delete column Col from B
	for (int i = 0; i < R - 1; ++i)
	{
		B[i].erase(B[i].begin() + Col);
	}
	/*DenseMatrix N = DenseMatrix::DenseMatrix(R - 1, C - 1);
	N.Matrix = B;
	N.Rows = R - 1;
	N.Cols = C - 1;
	PrintMatrix(&N); */
	return pow(-1, Row + Col) * DetHelper(B, R - 1, C - 1);
}

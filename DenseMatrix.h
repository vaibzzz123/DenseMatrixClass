#pragma once
#include <iostream>
#include <vector>
#include <cmath>

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
	static double Determinant(DenseMatrix *M);

	//inverse of a matrix
	bool IsInvertible(DenseMatrix *M);
	static DenseMatrix* MatrixInverse(DenseMatrix *M);

	//some IO functions
	static void PrintMatrix(DenseMatrix *M);
	static void PrintVector(std::vector<double>* V);

private:

};

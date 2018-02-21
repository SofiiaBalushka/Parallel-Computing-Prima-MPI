#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <time.h>
#include <mpi.h>
#include <vector>
#include <iostream>
using namespace std;

int ProcNum = 0; 
int ProcRank = 0; 


void RandomDataInitialization(double* pMatrix, int Size) {
	int i, j; 
	srand(unsigned(clock()));
	for (i = 0; i<Size; i++) {
		for (j = 0; j<Size; j++)
			pMatrix[i*Size + j] = rand() / double(1000);
	}
}


void ProcessInitialization(double* &pMatrix, double* &pResult, double* &pProcRows, double* &pProcResult, int &Size, int &RowNum) 
{
	int RestRows; 
	int i; 
	if (ProcRank == 0) {
		do {
			printf("\nEnter size of the initial objects: "); cout.flush();
			scanf("%d", &Size);
			if (Size < ProcNum) {
				printf("Size of the objects must be greater than number of processes!\n "); cout.flush();
			}
		} while (Size < ProcNum);
	}

	MPI_Bcast(&Size, 1, MPI_INT, 0, MPI_COMM_WORLD);
	RestRows = Size;
	for (i = 0; i<ProcRank; i++)
		RestRows = RestRows - RestRows / (ProcNum - i);
	RowNum = RestRows / (ProcNum - ProcRank);
	pResult = new double[Size*Size];
	pProcRows = new double[RowNum*Size];
	pProcResult = new double[RowNum*Size];
	if (ProcRank == 0) {
		pMatrix = new double[Size*Size];
		RandomDataInitialization(pMatrix, Size);
	}
}

void DataDistribution(double* pMatrix, double* pProcRows, int Size, int RowNum) 
{
	int *pSendNum; // the number of elements sent to the process
	int *pSendInd; // the index of the first data element sent to the process
	int RestColumns = Size; 
	
	pSendInd = new int[ProcNum];
	pSendNum = new int[ProcNum];
	
	// Define the disposition of the matrix rows for current process
	RowNum = (Size / ProcNum);
	pSendNum[0] = RowNum*Size;
	pSendInd[0] = 0;
	for (int i = 1; i<ProcNum; i++) {
		RestColumns -= RowNum;
		RowNum = RestColumns / (ProcNum - i);
		pSendNum[i] = RowNum*Size;
		pSendInd[i] = pSendInd[i - 1] + pSendNum[i - 1];
	}
	// Scatter the rows
	MPI_Scatterv(pMatrix, pSendNum, pSendInd, MPI_DOUBLE, pProcRows, pSendNum[ProcRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	
	delete[] pSendNum;
	delete[] pSendInd;
}


void ResultReplication(double* pProcResult, double* pResult, int Size,int RowNum) 
{
	int i; 
	int *pReceiveNum; // Number of elements, that current process sends
	int *pReceiveInd; // Index of the first element 
	int RestRows = Size; // Number of rows, that haven’t been distributed yet
	
	pReceiveNum = new int[ProcNum];
	pReceiveInd = new int[ProcNum];
	
	pReceiveInd[0] = 0;
	pReceiveNum[0] = Size / ProcNum;
	for (i = 1; i<ProcNum; i++) {
		RestRows -= pReceiveNum[i - 1];
		pReceiveNum[i] = RestRows / (ProcNum - i);
		pReceiveInd[i] = pReceiveInd[i - 1] + pReceiveNum[i - 1];
	}
	
		MPI_Allgatherv(pProcResult, pReceiveNum[ProcRank], MPI_DOUBLE, pResult, pReceiveNum, pReceiveInd, MPI_DOUBLE, MPI_COMM_WORLD);

	delete[] pReceiveNum;
	delete[] pReceiveInd;
}

void SerialResultCalculation(double* pMatrix, double* pResult, int Size) 
{
	vector<int> k (Size);
	for (int t = 0; t < Size; t++)
	{
		k[t] = 0;
	}
	int min = 0;
	for (int i = 0; i < Size; i++)
	{
		pResult[i] = 0;
		min = 10000000000;
		for (int j = 0; j < Size; j++)
		{
			if (pMatrix[i*Size + j] < min && pMatrix[i*Size + j] != 0)
			{
				k[j] = pMatrix[i*Size + j];
				min = pMatrix[i*Size + j];
			}
		}
	}
	for (int i = 0; i < Size; i++)
		for (int j = 0; j < Size; j++)
		{
			if (pMatrix[i*Size + j] != k[j])
				pResult[i*Size + j] = 0;
			else pResult[i*Size + j] = pMatrix[i*Size + j];
		}
}
void ParallelResultCalculation(double* pProcRows, double* pProcResult, int Size, int RowNum) 
{
	vector<int> k (Size);
	for (int t = 0; t < Size; t++)
	{
		k[t] = 0;
	}
	int min = 0;
	for (int i = 0; i < RowNum; i++)
	{
		pProcResult[i] = 0;
		min = 10000000000;
		for (int j = 0; j < Size; j++)
		{
			if (pProcRows[i*Size+j] < min && pProcRows[i*Size+j] != 0)
			{
				k[j] = pProcRows[i*Size + j];
				min = pProcRows[i*Size + j];
			}
		}
	}
	for (int i = 0; i < RowNum; i++)
		for (int j = 0; j < Size; j++)
		{
			if (pProcRows[i*Size + j] != k[j])
				pProcResult[i*Size + j] = 0;
			else pProcResult[i*Size + j] = pProcRows[i*Size + j];
		}
		
}


void PrintMatrix(double* pMatrix, int ColumnCount, int RowCount) 
{
	int i, j; // Loop variables
	for (i = 0; i<ColumnCount; i++) {
		for (j = 0; j<RowCount; j++)
			printf("%7.4f ", pMatrix[i*ColumnCount + j]);
		printf("\n");
	}
}

/*void PrintVector(double* pVector, int Size) 
{
	int i;
	for (i = 0; i<Size; i++)
		printf("%7.4f ", pVector[i]);
}*/

void TestDistribution(double* pMatrix, double* pProcColumns, int Size, int ColumnNum) 
{
	if (ProcRank == 0) {
		printf("Initial Matrix: \n");
		PrintMatrix(pMatrix, Size, Size);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	for (int i = 0; i<ProcNum; i++) {
		if (ProcRank == i) {
			printf("\nProcRank = %d \n", ProcRank);
			printf(" Matrix Stripe:\n");
			PrintMatrix(pProcColumns, ColumnNum, Size);
			
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
}

/*void TestPartialResults(double* pProcResult, int RowNum) {
	int i; // Loop variables
	for (i = 0; i<ProcNum; i++) {
		if (ProcRank == i) {
			printf("\nProcRank = %d \n Part of result vector: \n", ProcRank);
			//PrintVector(pProcResult, RowNum);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
}
*/
void TestResult(double* pMatrix, double* pResult, int Size) 
{
	
	double* pSerialResult;

	int equal = 0;
	int i; 
	if (ProcRank == 0) {
		pSerialResult = new double[Size*Size];
		SerialResultCalculation(pMatrix, pSerialResult, Size);
		for (i = 0; i<Size*Size; i++) {
			if (pResult[i] != pSerialResult[i])
				equal = 1;
		}
		if (equal == 1)
		{
			printf("The results of serial and parallel algorithms are NOT identical. Check your code."); cout.flush();
		}
		else
			printf("The results of serial and parallel algorithms are identical."); cout.flush();
	}
}


void ProcessTermination(double* pMatrix, double* pResult, double* pProcColumns, double* pProcResult) {
	if (ProcRank == 0)
		delete[] pMatrix;
	
	delete[] pResult;
	delete[] pProcColumns;
	delete[] pProcResult;
}

void main(int argc, char* argv[]) {
	double* pMatrix; // The first argument - initial matrix
	double* pResult; // Result matrix
	int Size; // Sizes of initial matrix 
	double* pProcRows;
	double* pProcResult;
	int RowNum;
	double Start, Finish, Duration;
	
	MPI_Init(&argc, &argv);
	
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	ProcessInitialization(pMatrix, pResult, pProcRows, pProcResult, Size, RowNum);
	Start = MPI_Wtime();
	DataDistribution(pMatrix, pProcRows,  Size, RowNum); 
	ParallelResultCalculation(pProcRows, pProcResult, Size, RowNum);
	ResultReplication(pProcResult, pResult, Size, RowNum);
	Finish = MPI_Wtime();
	Duration = Finish - Start;
	TestResult(pMatrix, pResult, Size);
	if (ProcRank == 0) {
		printf("Time of execution = %f\n", Duration); cout.flush();
	}
	ProcessTermination(pMatrix, pResult, pProcRows, pProcResult);
	MPI_Finalize();
}
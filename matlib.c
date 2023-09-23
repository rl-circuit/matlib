/*
    This library gives you the ability to manipulate matrices and their properties.
    Created by: Rares Luchian (lr.professional@outlook.com) in Anno Domini MMXXIII.
    This is FOSS. Enjoy! :))
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "approximations.h"

#define TRUE 1
#define FALSE 0

#define MEMERR "ERROR: there's not enough memory.\n"
#define NO_MTRX_FREE "There's no matrix to free.\n"
#define NO_MTRX_PRINT "There's no matrix to print.\n"
#define NO_MTRX_POPUL "There's no matrix to populate.\n"
#define NO_MTRX_CPY "There's no matrix to copy.\n"
#define MTRX_NOT_INV "Matrix not invertible.\n"
#define MTRX_EMPTY "ERROR: one or more matrices empty...\n"
#define MTRX_SPC_WRONG "ERROR: matrices spaces don't match.\n"
#define COLERROR "ERROR: index_1 or index_2 selected in switchCols() is greater than cols itself.\n"
#define ROWERROR "ERROR: index_1 or index_2 selected in switchRows() is greater than rows itself.\n"
#define ROW_ECH_ERROR "ERROR: Impossible to calculate row-echelon form. One column empty. The given matrix is not m*n!\n"

double** mallocMatrix(int rows, int cols)
{
    /*
        ACHTUNG! If malloc() fails in any circumstance and returns NULL,
        the whole subroutine mallocMatrix returns NULL like a regular malloc().
        It's up to you how to handle the NULL pointer in main().
    */
    
    size_t i, j;
    double** r = malloc(rows * sizeof(double*));
    
    if(r != NULL) {
        for(i = 0; i < rows; i++) {
            *(r + i) = malloc(cols * sizeof(double));
                        
            if(*(r + i) == NULL) {
                fprintf(stdout, MEMERR);
                for(j = 0; j < i; j++) {
                    free(*(r + j));
                }
                free(r);
                return NULL;
            }
        }
        return r;
    }
    else {
        fprintf(stdout, MEMERR);
        return NULL;
    }
}

void freeMatrix(double** r, int rows)
{
    size_t i;
    
    if(r != NULL) {
        for(i = 0; i < rows; i++) {
            free(*(r + i));
        }
        free(r);
    }
    else {
        fprintf(stdout, NO_MTRX_FREE);
        return;
    }
}

void printMatrix(double** r, int rows, int cols)
{
    size_t i, j;
    
    if(r != NULL) {
        for(i = 0; i < rows; i++) {
            for(j = 0; j < cols; j++) {
                fprintf(stdout, "%8.5f   ", r[i][j]);
            }
            fprintf(stdout, "\n");
        }
    }
    else {
        fprintf(stdout, NO_MTRX_PRINT);
        return;
    }
}

void populateMatrix(double** r, int rows, int cols)
{
    size_t i, j;
    
    if(r != NULL) {
        for(i = 0; i < rows; i++) {
            for(j = 0; j < cols; j++) {
                fprintf(stdout, "Type element [%ld][%ld] and press ENTER: ", i+1, j+1);
                fscanf(stdin, "%lf", &r[i][j]);
            }
        }
    }
    else {
        fprintf(stdout, NO_MTRX_POPUL);
        return;
    }
}

double** copyMatrix(double** r, int rows, int cols)
{
    /*
        ACHTUNG! If there's no matrix to copy or if mallocMatrix() fails to allocate,
        the subroutine returns NULL. It's up to you how to handle the null pointer
        in main().
    */
    
    double **copy;
    size_t i, j;
    
    if(r != NULL) {
        copy = mallocMatrix(rows, cols);
        if(copy != NULL) {
            for(i = 0; i < rows; i++) {
                for(j = 0; j < cols; j++) {
                    copy[i][j] = r[i][j];
                }
            }
            return copy;
        }
        else {
            fprintf(stdout, "ERROR: mallocMatrix() failed in copyMatrix().\n");
            return NULL;
        }
    }
    else {
        fprintf(stdout, NO_MTRX_CPY);
        return NULL;
    }
}

double** createIdentity(int rows)
{
    size_t i, j;
    
    double** I = mallocMatrix(rows, rows);
    if(I != NULL) {
        for(i = 0; i < rows; i++) {
            for(j = 0; j < rows; j++) {
                if(i == j) I[i][j] = 1;
                else {
                    I[i][j] = 0;
                }
            }
        }
        return I;
    }
    else {
        fprintf(stdout, "ERROR: mallocMatrix() failed in createIdentity().\n");
        return NULL;        
    }
}

void ifSwitchableRows(double** r, int rows, int index_1, int index_2, int* value)
{
    /*
        ACHTUNG! index_1 refers to which row the zero is found. Index_2 refers to
        which column the zero is found. Value is the possible row index to switch with.
        If no candidates are found (the entire row is made of 0), then value = -1.
    */
    
    int i, count = 0;

    for(i = index_1; i < rows; i++) {
        approxZero(&r[i][index_2]);
        if(r[i][index_2] == 0) count++;
        else {
            *value = i;
            break;
        }
    }
    if(count == (rows - index_1)) *value = -1;
}

void ifSwitchableRowsReverse(double** r, int rows, int index_1, int index_2, int* value)
{
    /*
        ACHTUNG! index_1 refers to which row the zero is found. Index_2 refers to
        which column the zero is found. Value is the possible row index to switch with.
        If no candidates are found (the entire row is made of 0), then value = -1.
    */
    
    int i, count = 0;

    for(i = index_1; i >= 0; i--) {
        approxZero(&r[i][index_2]);
        if(r[i][index_2] == 0) count++;
        else {
            *value = i;
            break;
        }
    }
    if(count == index_1) *value = -1;
}

void ifSwitchableCols(double** r, int cols, int index_1, int index_2, int* value)
{
    /*
        ACHTUNG! index_1 refers to which row the zero is found. Index_2 refers to
        which column the zero is found. Value is the possible column index to switch with.
        If no candidates are found (the entire column is made of 0), then value = -1.
    */
    
    int i, count = 0;
    
    for(i = index_2; i < cols; i++) {
        approxZero(&r[index_1][i]);
        if(r[index_1][i] == 0) count++;
        else {
            *value = i;
            break;
        }
    }
    if(count == (cols - index_2)) *value = -1;
}

void switchCols(double** r, int rows, int cols, int index_1, int index_2)
{
    size_t i;
    double temp;
    
    if(index_1 > cols || index_2 > cols) {
        fprintf(stdout, COLERROR);
        return;
    }
    else {
        if(r) {
            for(i = 0; i < rows; i++) {
                temp = r[i][index_1];
                r[i][index_1] = r[i][index_2];
                r[i][index_2] = temp;
            }
        }
        else {
            fprintf(stdout, MTRX_EMPTY);
            return;
        }
    }
}

void switchRows(double** r, int rows, int cols, int index_1, int index_2)
{
    size_t i;
    double temp;
    
    if(index_1 > rows || index_2 > rows) {
        fprintf(stdout, ROWERROR);
        return;
    }
    else {
        if(r) {
            for(i = 0; i < cols; i++) {
                temp = r[index_1][i];
                r[index_1][i] = r[index_2][i];
                r[index_2][i] = temp;
            }
        }
        else {
            fprintf(stdout, MTRX_EMPTY);
            return;
        }
    }
}

void gaussianUpper(double** r, int rows, int cols, int* switchedRows)
{
    size_t i, j, k;
    double temp, factor;
    int index = -2;
        
    if(r == NULL) {
        fprintf(stdout, MTRX_EMPTY);
        return;
    }
    else {
        for(k = 0; k < rows-1; k++) {
            approxZero(&r[k][k]);
            temp = r[k][k];
            
            if(temp == 0) {
                ifSwitchableRows(r, rows, k, k, &index);
                if(index == -1) {
                    fprintf(stdout, ROW_ECH_ERROR);
                    return;
                }
                else {
                    switchRows(r, rows, cols, k, index);
                    (*switchedRows)++;
                    temp = r[k][k];
                }
            }
            
            for(i = 0; i < rows-k-1; i++) {
                factor = -r[k+i+1][k] / temp;
                for(j = 0; j < cols; j++) {
                    r[k+i+1][j] += factor * r[k][j];
                    approxZero(&r[k+i+1][j]);
                }
            }
        }
    }    
}

void gaussianLower(double** r, int rows, int cols, int* switchedRows)
{
    size_t i, j, k;
    double temp, factor;
    int index = -2;    
    
    if(r == NULL) {
        fprintf(stdout, MTRX_EMPTY);
        return;
    }
    else {
        for(k = 0; k < rows-1; k++) {
            approxZero(&r[rows-k-1][cols-k-1]);
            temp = r[rows-k-1][cols-k-1];

            if(temp == 0) {
                ifSwitchableRowsReverse(r, rows, rows-k-1, cols-k-1, &index);
                if(index == -1) {
                    fprintf(stdout, ROW_ECH_ERROR);
                    return;
                }
                else {
                    switchRows(r, rows, cols, rows-k-1, index);
                    (*switchedRows)++;
                    temp = r[rows-k-1][cols-k-1];
                }
            }
            
            for(i = 0; i < rows-k-1; i++) {
                factor = -r[rows-k-i-2][cols-k-1] / temp;
                for(j = 0; j < cols; j++) {
                    r[rows-k-i-2][cols-j-1] += factor * r[rows-k-1][cols-j-1];
                    approxZero(&r[rows-k-i-2][cols-j-1]);
                }
            }
        }
    }
}

void normalize(double** r, int rows)
{
    /*
        ACHTUNG! Only for square matrices!
    */
    
    size_t k, j;
    double temp;
    
    if(r) {
        for(k = 0; k < rows; k++) {
            temp = r[k][k];
            for(j = 0; j < rows; j++) {
                r[k][j] /= temp;
            }
        }
    }
    else {
        fprintf(stdout, MTRX_EMPTY);
        return;
    }
}

double** transpose(double** r, int rows, int cols, int* tRows, int* tCols)
{
    /*
        ACHTUNG! Remember to free transposed in main() after use!
    */
    
    size_t i, j;
    
    if(r == NULL) {
        fprintf(stdout, MTRX_EMPTY);
        return NULL;
    }
    else {
        double** transposed = mallocMatrix(cols, rows);
        
        if(transposed != NULL) {
            *tRows = cols;
            *tCols = rows;
            
            for(i = 0; i < rows; i++) {
                for(j = 0; j < cols; j++) {
                    transposed[j][i] = r[i][j];
                }
            }
            return transposed;
        }
        else {
            fprintf(stdout, "ERROR: mallocMatrix() failed in transpose().\n");
            return NULL;
        }
    }
}


int rank(double** r, int rows, int cols)
{
    int count = 0, switched = 0, min, i = 0;
    double** copy = copyMatrix(r, rows, cols);
    
    if(r != NULL && copy != NULL) {
        gaussianUpper(copy, rows, cols, &switched);
        
        if(rows < cols) min = rows;
        else min = cols;
        while(i < min) {
            if(copy[i][i] != 0) count++;
            i++;
        }
        freeMatrix(copy, rows);
        return count;
    }
    else {
        fprintf(stdout, MTRX_EMPTY);
        return -1;
    }
}

double determinant(double** r, int rows)
{
    int switches = 0, i;
    double det = 1.0, **copy; 
    copy = copyMatrix(r, rows, rows);
    
    if(copy) {
        gaussianUpper(copy, rows, rows, &switches);
        for (i = 0; i < rows; i++) {
            det *= copy[i][i];
        }
        freeMatrix(copy, rows);
        return pow(-1, switches) * det;
    }
    else {
        fprintf(stdout, "ERROR: copyMatrix() failed in determinant(). The determinant is incorrect.\n");
        return -1;
    }
}

double** matrixInversion(double** r, int rows)
{
    /*
        ACHTUNG! Remember to free I in main() after use!
    */
    
    int i, j, k;
    double temp, factor, **I, **A;
    int index = -2, rango, switchedRows = 0;
    
    if(r) {
        rango = rank(r, rows, rows);
        if(rango == rows) {
            I = createIdentity(rows);
            A = copyMatrix(r, rows, rows);
            
            if(I) {
                if(A) {
                    for(k = 0; k < rows-1; k++) {
                        approxZero(&A[k][k]);
                        temp = A[k][k];
                    
                        if(temp == 0) {
                            ifSwitchableRows(A, rows, k, k, &index);
                            if(index == -1) {
                                fprintf(stdout, ROW_ECH_ERROR);
                                return NULL;
                            }
                            else {
                                switchRows(A, rows, rows, k, index);
                                switchRows(I, rows, rows, k, index);
                                switchedRows++;
                                temp = A[k][k];
                            }
                        }
                    
                        for(i = 0; i < rows-k-1; i++) {
                            factor = -A[k+i+1][k] / temp;
                            for(j = 0; j < rows; j++) {
                                A[k+i+1][j] += factor * A[k][j];
                                approxZero(&A[k+i+1][j]);
                                I[k+i+1][j] += factor * I[k][j];
                                approxZero(&I[k+i+1][j]);
                            }
                        }
                    }
                    
                    for(k = 0; k < rows-1; k++) {
                        approxZero(&A[rows-k-1][rows-k-1]);
                        temp = A[rows-k-1][rows-k-1];

                        if(temp == 0) {
                            ifSwitchableRowsReverse(A, rows, rows-k-1, rows-k-1, &index);
                            if(index == -1) {
                                fprintf(stdout, ROW_ECH_ERROR);
                                return NULL;
                            }
                            else {
                                switchRows(A, rows, rows, rows-k-1, index);
                                switchRows(I, rows, rows, rows-k-1, index);
                                switchedRows++;
                                temp = A[rows-k-1][rows-k-1];
                            }
                        }
            
                        for(i = 0; i < rows-k-1; i++) {
                            factor = -A[rows-k-i-2][rows-k-1] / temp;
                            for(j = 0; j < rows; j++) {
                                A[rows-k-i-2][rows-j-1] += factor * A[rows-k-1][rows-j-1];
                                approxZero(&A[rows-k-i-2][rows-j-1]);
                                I[rows-k-i-2][rows-j-1] += factor * I[rows-k-1][rows-j-1];
                                approxZero(&I[rows-k-i-2][rows-j-1]);
                            }
                        }
                    }
                    
                    for(k = 0; k < rows; k++) {
                        temp = A[k][k];
                        for(j = 0; j < rows; j++) {
                            A[k][j] /= temp;
                            I[k][j] /= temp;
                        }
                    }
                    return I;
                }
                else {
                    fprintf(stdout, MEMERR);
                    return NULL;
                }
            }
            else {
                fprintf(stdout, "ERROR: createIdentity() failed in matrixInversion()\n");
                return NULL;
            }
        }
        else {
            fprintf(stdout, MTRX_NOT_INV);
            return NULL;
        }
    }
    else {
        fprintf(stdout, MTRX_EMPTY);
        return NULL;
    }
}

double** multiply(double** A, double** B, int rowsA, int colsA, int rowsB, int colsB, int* rowsProd, int* colsProd)
{
	if(A == NULL || B == NULL) {
        fprintf(stdout, MTRX_EMPTY);
        return NULL;
	}
	else if(colsA != rowsB) {
	    fprintf(stdout, MTRX_SPC_WRONG);
        return NULL;
	}
	else {
	    size_t i, j, k;
	    double** product = mallocMatrix(rowsA, colsB);
	    
	    if(product != NULL) {
	        *rowsProd = rowsA;
	        *colsProd = colsB;
	        
	        for(i = 0; i < rowsA; i++) {
	            for(j = 0; j < colsB; j++) {
	                product[i][j] = 0;
	                for(k = 0; k < colsA; k++) {
	                    product[i][j] += A[i][k] * B[k][j];
	                    approxZero(&product[i][j]);
	                }
	            }
	        }
	        return product;
	        /* ACHTUNG! Remember to free product in main() after use! */
	    }
	    else {
            fprintf(stdout, "ERROR: mallocMatrix() failed in product().\n");
            return NULL;
	    }
	}
}

double** add(double** A, double** B, int rowsA, int colsA, int rowsB, int colsB)
{
    if(A == NULL || B == NULL) {
        fprintf(stdout, MTRX_EMPTY);
        return NULL;
    }
    else if(rowsA != rowsB || colsA != colsB) {
        fprintf(stdout, MTRX_SPC_WRONG);
        return NULL;
    }
    else {
        size_t i, j;
        double** result = mallocMatrix(rowsA, colsA);
        
        if(result != NULL) {
            for(i = 0; i < rowsA; i++) {
                for(j = 0; j < colsA; j++) {
                    result[i][j] = A[i][j] + B[i][j];
                    approxZero(&result[i][j]);
                }
            }
	        return result;
	        /* ACHTUNG! Remember to free result in main() after use! */
        }
        else {
            fprintf(stdout, "ERROR: mallocMatrix() failed in add().\n");
            return NULL;
        }
    }
}

void scalarMultiply(double** A, int rowsA, int colsA, double scalar)
{
    if(A == NULL) {
        fprintf(stdout, MTRX_EMPTY);
        return;
    }
    else {
        size_t i, j;        
        for(i = 0; i < rowsA; i++) {
            for(j = 0; j < colsA; j++) {
                *(*(A + i) + j) *= scalar;
                approxZero(&A[i][j]);
            }
        }
    }
}

#ifndef MATLIB_H

    #define MATLIB_H
    
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
    
    double** mallocMatrix(int rows, int cols);
    void freeMatrix(double** r, int rows);
    void printMatrix(double** r, int rows, int cols);
    void populateMatrix(double** r, int rows, int cols);
    double** copyMatrix(double** r, int rows, int cols);
    double** createIdentity(int rows);
    void ifSwitchableRows(double** r, int rows, int index_1, int index_2, int* value);
    void ifSwitchableRowsReverse(double** r, int rows, int index_1, int index_2, int* value);
    void ifSwitchableCols(double** r, int cols, int index_1, int index_2, int* value);
    void switchCols(double** r, int rows, int cols, int index_1, int index_2);
    void switchRows(double** r, int rows, int cols, int index_1, int index_2);
    void gaussianUpper(double** r, int rows, int cols, int* switchedRows);
    void gaussianLower(double** r, int rows, int cols, int* switchedRows);
    void normalize(double** r, int rows);
    double** transpose(double** r, int rows, int cols, int* tRows, int* tCols);
    int rank(double** r, int rows, int cols);
    double determinant(double** r, int rows);
    double** matrixInversion(double** r, int rows);
    double** multiply(double** A, double** B, int rowsA, int colsA, int rowsB, int colsB, int* rowsProd, int* colsProd);
    double** add(double** A, double** B, int rowsA, int colsA, int rowsB, int colsB);
    void scalarMultiply(double** A, int rowsA, int colsA, double scalar);
    
#endif

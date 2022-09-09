#include <stdio.h>
#include <math.h>

// Struct representing a 3D vector
typedef struct {
    int x;
    int y;
    int z;
} Vector;

// Struct representing a 3D matrix
typedef struct {
    Vector c1;
    Vector c2;
    Vector c3;
} Matrix;

// Vector operations
double getVectorMagnitude(Vector*);
int areVectorsEqual(Vector*, Vector*);
void scaleVector(Vector*, int);
void addVectors(Vector*, Vector*, Vector*);
double dotProduct(Vector*, Vector*);
void crossProduct(Vector*, Vector*, Vector*);

// Matrix operations
void multiplyMatrixTimesVector(Matrix*, Vector*, Vector*);
int areMatricesEqual(Matrix*, Matrix*);
void addMatrices(Matrix*, Matrix*, Matrix*);
void multiplyMatrices(Matrix*, Matrix*, Matrix*);
void transposeMatrix(Matrix*, Matrix*);
double getDeterminant(Matrix*);
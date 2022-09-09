#include "linalg.h"

// Returns the length of the given vector
double getVectorMagnitude(Vector* vector) {
    return sqrt(
        pow(vector->x, 2) + 
        pow(vector->y,  2) +
        pow(vector->z, 2)
    );
}

// Returns 1 if given vectors are equal, 0 otherwise
int areVectorsEqual(Vector* vector1, Vector* vector2) {
    return vector1->x == vector2->x && vector1->y == vector2->y && vector1->z == vector2->z;
}

// Multiplies the given vector by the given scalar
void scaleVector(Vector* vector, int scalar) {
    vector->x *= scalar;
    vector->y *= scalar;
    vector->z *= scalar;
}

// Adds the two vectors and places the sum in the sum vector
void addVectors(Vector* vector1, Vector* vector2, Vector* sum) {
    sum->x = vector1->x + vector2->x;
    sum->y = vector1->y + vector2->y;
    sum->z = vector1->z + vector2->z;
}

// Returns the dot product of the given vectors
double dotProduct(Vector* vector1, Vector* vector2) {
    return vector1->x * vector2->x +
           vector1->y * vector2->y +
           vector1->z * vector2->z;
}

// Calculates the cross product of vector1 and vector2 and puts the result in the
// product vector
void crossProduct(Vector* vector1, Vector* vector2, Vector* product) {
    product->x = vector1->y * vector2->z - vector1->z * vector2->y;
    product->y = -(vector1->x * vector2->z - vector1->z * vector2->x);
    product->z = vector1->x * vector2->y - vector1->y * vector2->x;
}

// Multiples the matrix by the vector and puts the result in the product vector
void multiplyMatrixTimesVector(Matrix* matrix, Vector* vector, Vector* product) {
    // Vectors to be scaled by elements of given vector
    Vector scaledVector1 = matrix->c1;
    Vector scaledVector2 = matrix->c2;
    Vector scaledVector3 = matrix->c3;
    
    scaleVector(&scaledVector1, vector->x);
    scaleVector(&scaledVector2, vector->y);
    scaleVector(&scaledVector3, vector->z);

    // Calculating result vector after linear transformation
    product->x = scaledVector1.x + scaledVector2.x + scaledVector3.x;
    product->y = scaledVector1.y + scaledVector2.y + scaledVector3.y;
    product->z = scaledVector1.z + scaledVector2.z + scaledVector3.z;
}

// Returns 1 if given matrices are equal, 0 otherwise
int areMatricesEqual(Matrix* matrix1, Matrix* matrix2) {
    return areVectorsEqual(&matrix1->c1, &matrix2->c1) && areVectorsEqual(&matrix1->c2, &matrix2->c2) && areVectorsEqual(&matrix1->c3, &matrix2->c3);
}

// Adds the two matrices together and puts the result into sum
void addMatrices(Matrix* matrix1, Matrix* matrix2, Matrix* sum) {
    addVectors(&matrix1->c1, &matrix2->c1, &sum->c1);
    addVectors(&matrix1->c2, &matrix2->c2, &sum->c2);
    addVectors(&matrix1->c3, &matrix2->c3, &sum->c3);
}

// Multiplies to the two matrices together and puts result into product
void multiplyMatrices(Matrix* matrix1, Matrix* matrix2, Matrix* product) {
    multiplyMatrixTimesVector(matrix1, &matrix2->c1, &product->c1);
    multiplyMatrixTimesVector(matrix1, &matrix2->c2, &product->c2);
    multiplyMatrixTimesVector(matrix1, &matrix2->c3, &product->c3);
}

// Transposes the given matrix and puts the result into transpose
void transposeMatrix(Matrix* matrix, Matrix* transpose) {
    transpose->c1 = (Vector){ matrix->c1.x, matrix->c2.x, matrix->c3.x };
    transpose->c2 = (Vector){ matrix->c1.y, matrix->c2.y, matrix->c3.y };
    transpose->c3 = (Vector){ matrix->c1.z, matrix->c2.z, matrix->c3.z };
}

// Returns the determinant of the given matrix
double getDeterminant(Matrix* matrix) {
    // Calculate each component of the determinate individually
    double iComponent = matrix->c1.x * (matrix->c2.y * matrix->c3.z - matrix->c3.y * matrix->c2.z);
    double jComponent = matrix->c2.x * (matrix->c1.y * matrix->c3.z - matrix->c3.y * matrix->c1.z);
    double kComponent = matrix->c3.x * (matrix->c1.y * matrix->c2.z - matrix->c2.y * matrix->c1.z);
    
    return iComponent - jComponent + kComponent;
}
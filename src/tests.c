#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "linalg.h"

#define FLOAT_THRESH 0.00000001 // Threshold for how close floats need to be to assert them equal

void testVectorMagnitude() {
    assert(fabs(getVectorMagnitude(&(Vector){0, 0, 0})) == 0);
    assert(fabs(getVectorMagnitude(&(Vector){4, 1, 1}) - 3 * sqrt(2)) < FLOAT_THRESH);
    assert(fabs(getVectorMagnitude(&(Vector){3, 4, 5}) - 5 * sqrt(2)) < FLOAT_THRESH);

    // Negative tests
    assert(fabs(getVectorMagnitude(&(Vector){-4, 1, 1}) - 3 * sqrt(2)) < FLOAT_THRESH);
    assert(fabs(getVectorMagnitude(&(Vector){-3, -4, 5}) - 5 * sqrt(2)) < FLOAT_THRESH);
}

void testAreVectorsEqual() {
    assert(!areVectorsEqual(&(Vector){0, 0, 0}, &(Vector){1, 5, 6}));
    assert(areVectorsEqual(&(Vector){2, 1, 8}, &(Vector){2, 1, 8}));
    assert(areVectorsEqual(&(Vector){0, 0, 0}, &(Vector){0, 0, 0}));
    assert(areVectorsEqual(&(Vector){-5, 1, 2}, &(Vector){-5, 1, 2}));

    assert(!areVectorsEqual(&(Vector){5, 1, 2}, &(Vector){5, 1, 1}));
    assert(!areVectorsEqual(&(Vector){5, 1, 2}, &(Vector){5, 0, 2}));
    assert(!areVectorsEqual(&(Vector){5, 1, 2}, &(Vector){0, 1, 2}));
}

void testScaleVector() {
    Vector vector1;
    scaleVector(&(Vector){1, 2, 4}, 0, &vector1);
    assert(areVectorsEqual(&vector1, &(Vector){0, 0, 0}));

    Vector vector2;
    scaleVector(&(Vector){0, 0, 0}, 5, &vector2);
    assert(areVectorsEqual(&vector2, &(Vector){0, 0, 0}));

    Vector vector3;
    scaleVector(&(Vector){5, 4, 7}, 1, &vector3);
    assert(areVectorsEqual(&vector3, &(Vector){5, 4, 7}));

    Vector vector4;
    scaleVector(&(Vector){1, 1, 1}, 2, &vector4);
    assert(areVectorsEqual(&vector4, &(Vector){2, 2, 2}));

    Vector vector5;
    scaleVector(&(Vector){4, 1, 2}, 3, &vector5);
    assert(areVectorsEqual(&vector5, &(Vector){12, 3, 6}));

    Vector vector6;
    scaleVector(&(Vector){6, 9, 1}, -2, &vector6);
    assert(areVectorsEqual(&vector6, &(Vector){-12, -18, -2}));
}

void testAddVectors() {
    Vector sum1;
    addVectors(&(Vector){0, 0, 0}, &(Vector){1, 5, 3}, &sum1);
    assert(areVectorsEqual(&sum1, &(Vector){1, 5, 3}));

    Vector sum2;
    addVectors(&(Vector){1, 3, 5}, &(Vector){2, 2, 2}, &sum2);
    assert(areVectorsEqual(&sum2, &(Vector){3, 5, 7}));

    Vector sum3;
    addVectors(&(Vector){0, 0, 0}, &(Vector){0, 0, 1}, &sum3);
    assert(areVectorsEqual(&sum3, &(Vector){0, 0, 1}));
    
    Vector sum4;
    addVectors(&(Vector){0, 0, 0}, &(Vector){0, 1, 0}, &sum4);
    assert(areVectorsEqual(&sum4, &(Vector){0, 1, 0}));

    Vector sum5;
    addVectors(&(Vector){0, 0, 0}, &(Vector){1, 0, 0}, &sum5);
    assert(areVectorsEqual(&sum5, &(Vector){1, 0, 0}));
}

void testDotProduct() {
    assert(fabs(dotProduct(&(Vector){0, 0, 0}, &(Vector){0, 0, 0})) < FLOAT_THRESH);
    assert(fabs(dotProduct(&(Vector){1, 5, 6}, &(Vector){4, 2, 1}) - 20) < FLOAT_THRESH);
    assert(fabs(dotProduct(&(Vector){1, 5, 6}, &(Vector){0, 0, 0}) - 0) < FLOAT_THRESH);
    assert(fabs(dotProduct(&(Vector){7, 0, 0}, &(Vector){0, 7, 0}) - 0) < FLOAT_THRESH);
    assert(fabs(dotProduct(&(Vector){7, -4, 0}, &(Vector){2, 7, -1}) + 14) < FLOAT_THRESH);
    assert(fabs(dotProduct(&(Vector){0, -1, 3}, &(Vector){5, 5, -5}) + 20) < FLOAT_THRESH);
}

void testCrossProduct() {
    Vector product1;
    crossProduct(&(Vector){0, 0, 0}, &(Vector){0, 0, 0}, &product1);
    assert(areVectorsEqual(&product1, &(Vector){0, 0, 0}));

    Vector product2;
    crossProduct(&(Vector){4, -1, 2}, &(Vector){5, 5, 2}, &product2);
    assert(areVectorsEqual(&product2, &(Vector){-12, 2, 25}));

    Vector product3;
    crossProduct(&(Vector){-1, 1, 1}, &(Vector){1, -1, 1}, &product3);
    assert(areVectorsEqual(&product3, &(Vector){2, 2, 0}));
}

void testMatricesEqual() {
    assert(areMatricesEqual(
        &(Matrix){
            (Vector){1, 4, 7},
            (Vector){2, 5, 8},
            (Vector){3, 6, 9}
        },
        &(Matrix){
            (Vector){1, 4, 7},
            (Vector){2, 5, 8},
            (Vector){3, 6, 9}
        }
    ));

    assert(!areMatricesEqual(
        &(Matrix){
            (Vector){1, 4, 7},
            (Vector){2, 5, 8},
            (Vector){3, 6, 9}
        },
        &(Matrix){
            (Vector){9, 6, 3},
            (Vector){8, 5, 2},
            (Vector){7, 4, 1}
        }
    ));

    assert(!areMatricesEqual(
        &(Matrix){
            (Vector){1, 4, 7},
            (Vector){2, 5, 8},
            (Vector){3, 6, 9}
        },
        &(Matrix){
            (Vector){0, 4, 7},
            (Vector){2, 5, 8},
            (Vector){3, 6, 9}
        }
    ));

    assert(!areMatricesEqual(
        &(Matrix){
            (Vector){1, 4, 7},
            (Vector){2, 5, 8},
            (Vector){3, 6, 9}
        },
        &(Matrix){
            (Vector){1, 4, 7},
            (Vector){3, 5, 8},
            (Vector){3, 6, 9}
        }
    ));

    assert(!areMatricesEqual(
        &(Matrix){
            (Vector){1, 4, 7},
            (Vector){2, 5, 8},
            (Vector){3, 6, 9}
        },
        &(Matrix){
            (Vector){1, 4, 7},
            (Vector){2, 5, 8},
            (Vector){0, 6, 9}
        }
    ));
}

void testMatrixTimesVector() {
    Vector product1;
    multiplyMatrixTimesVector(&(Matrix){
        (Vector){6, 2, 8},
        (Vector){3, -2, 5},
        (Vector){1, 7, 9},
    }, &(Vector){0, 0, 0}, &product1);
    assert(areVectorsEqual(&product1, &(Vector){0, 0, 0}));

    Vector product2;
    multiplyMatrixTimesVector(&(Matrix){
        (Vector){8, 2, 7},
        (Vector){1, 1, -8},
        (Vector){1, 7, 3},
    }, &(Vector){1, 5, 9}, &product2);
    assert(areVectorsEqual(&product2, &(Vector){22, 70, -6}));

    Vector product3;
    multiplyMatrixTimesVector(&(Matrix){
        (Vector){-1, -17, 40},
        (Vector){10, 3, 12},
        (Vector){11, 35, 13},
    }, &(Vector){9, -4, -6}, &product3);
    assert(areVectorsEqual(&product3, &(Vector){-115, -375, 234}));
}

void testAddMatrices() {
    Matrix sum1;
    addMatrices(
        &(Matrix){
            (Vector){1, 4, 7},
            (Vector){2, 5, 8},
            (Vector){3, 6, 9}
        },
        &(Matrix){
            (Vector){9, 6, 3},
            (Vector){8, 5, 2},
            (Vector){7, 4, 1}
        },
        &sum1
    );
    assert(areMatricesEqual(&sum1,
        &(Matrix){
            (Vector){10, 10, 10},
            (Vector){10, 10, 10},
            (Vector){10, 10, 10}
        }
    ));
}

void testMultiplyMatrices() {
    Matrix product1;
    multiplyMatrices(
    &(Matrix){
        (Vector){6, -1, 2},
        (Vector){9, 2, 4},
        (Vector){-4, 7, 3}
    },
    &(Matrix){
        (Vector){2, 8, 1},
        (Vector){12, -2, 3},
        (Vector){11, 7, 3}
    },
    &product1);

    assert(areMatricesEqual(&product1,
        &(Matrix){
            (Vector){80, 21, 39},
            (Vector){42, 5, 25},
            (Vector){117, 24, 59}
        }
    ));

    Matrix product2;
    multiplyMatrices(
    &(Matrix){
        (Vector){1, 2, 3},
        (Vector){8, 8, 6},
        (Vector){9, 2, 1}
    },
    &(Matrix){
        (Vector){-1, 0, 4},
        (Vector){8, 7, -10},
        (Vector){0, 4, 9}
    },
    &product2);

    assert(areMatricesEqual(&product2,
        &(Matrix){
            (Vector){35, 6, 1},
            (Vector){-26, 52, 56},
            (Vector){113, 50, 33}
        }
    ));

    Matrix product3;
    multiplyMatrices(
    &(Matrix){
        (Vector){1, 4, 7},
        (Vector){2, 5, 8},
        (Vector){3, 6, 9}
    },
    &(Matrix){
        (Vector){10, 13, 16},
        (Vector){11, 14, 17},
        (Vector){12, 15, 18}
    },
    &product3);

    assert(areMatricesEqual(&product3,
        &(Matrix){
            (Vector){84, 201, 318},
            (Vector){90, 216, 342},
            (Vector){96, 231, 366}
        }
    ));
}

void testTranspose() {
    Matrix transpose1;
    transposeMatrix(
        &(Matrix){
            (Vector){1, 2, 3},
            (Vector){4, 5, 6},
            (Vector){7, 8, 9}
        },
        &transpose1
    );

    assert(areMatricesEqual(&transpose1,
        &(Matrix){
            (Vector){1, 4, 7},
            (Vector){2, 5, 8},
            (Vector){3, 6, 9}
        }
    ));

    Matrix transpose2;
    transposeMatrix(
        &(Matrix){
            (Vector){1, 4, 7},
            (Vector){2, 5, 8},
            (Vector){3, 6, 9}
        },
        &transpose2
    );

    assert(areMatricesEqual(&transpose2,
        &(Matrix){
            (Vector){1, 2, 3},
            (Vector){4, 5, 6},
            (Vector){7, 8, 9}
        }
    ));
}

void testDeterminant() {
    assert(fabs(getDeterminant(
        &(Matrix){
            (Vector){-1, 2, 6},
            (Vector){2, 5, 8},
            (Vector){6, 7, 0}
        }
    ) - 56) < FLOAT_THRESH);

    assert(fabs(getDeterminant(
        &(Matrix){
            (Vector){1, 2, 3},
            (Vector){4, 5, 6},
            (Vector){7, 8, 9}
        }
    )) < FLOAT_THRESH);

    assert(fabs(getDeterminant(
        &(Matrix){
            (Vector){0, 0, 0},
            (Vector){0, 0, 0},
            (Vector){0, 0, 0},
        }
    )) < FLOAT_THRESH);

    assert(fabs(getDeterminant(
        &(Matrix){
            (Vector){1, 0, 0},
            (Vector){0, 1, 0},
            (Vector){0, 0, 1},
        }
    ) - 1) < FLOAT_THRESH);

    assert(fabs(getDeterminant(
        &(Matrix){
            (Vector){1, 2, 0},
            (Vector){2, 1, 2},
            (Vector){0, 0, 1},
        }
    ) + 3) < FLOAT_THRESH);
}

// Runs unit tests for the linalg functions
int main() {
    // Tests for vector operations
    testVectorMagnitude();
    testAreVectorsEqual();
    testScaleVector();
    testAddVectors();
    testDotProduct();
    testCrossProduct();

    // Tests for matrix operations
    testMatricesEqual();
    testMatrixTimesVector();
    testAddMatrices();
    testMultiplyMatrices();
    testTranspose();
    testDeterminant();

    printf("All tests ran successfully.\n");
    return 0;
}
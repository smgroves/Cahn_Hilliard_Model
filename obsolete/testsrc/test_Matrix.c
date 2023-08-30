#include <stdio.h>
#include <stdlib.h>
#include "unity.h" // Unity testing framework

// Include the header file for the function being tested
#include "C/CHsolver.h" // Replace with the actual header file name

// Define test cases
void test_dmatrix_allocation(void)
{
    long nrl = 1, nrh = 3, ncl = 1, nch = 4;
    double **matrix = dmatrix(nrl, nrh, ncl, nch);

    // Check if the allocation was successful
    TEST_ASSERT_NOT_NULL(matrix);

    // Check a few individual matrix elements
    TEST_ASSERT_EQUAL_DOUBLE(0.0, matrix[1][1]);
    TEST_ASSERT_EQUAL_DOUBLE(0.0, matrix[3][4]);

    // Free the allocated memory
    free(matrix[nrl]);
    free(matrix + nrl);
}

void test_dmatrix_invalid_indices(void)
{
    // Test with invalid indices
    long nrl = 3, nrh = 1, ncl = 4, nch = 1;
    double **matrix = dmatrix(nrl, nrh, ncl, nch);

    // Check if the function returns NULL for invalid indices
    TEST_ASSERT_NULL(matrix);
}

// Define the test runner function
int main(void)
{
    UNITY_BEGIN();

    // Run the test cases
    RUN_TEST(test_dmatrix_allocation);
    RUN_TEST(test_dmatrix_invalid_indices);

    return UNITY_END();
}
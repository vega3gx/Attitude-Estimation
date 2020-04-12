#include "MatrixMath.h"
#include <stdio.h>
#include <math.h>

void MatrixMultiply(double mat1[3][3], double mat2[3][3], double result[3][3]) {
    double dpod;
    int j, k;
    dpod = 0;
    for (j = 0; j < 3; j = j + 1) { //This is the same combo for loops I will use in several funtions
        for (k = 0; k < 3; k = k + 1) {
            dpod = (mat1[j][0] * mat2[0][k])+(mat1[j][1] * mat2[1][k])+(mat1[j][2] * mat2[2][k]);
            result[j][k] = dpod;
            dpod = 0;
        }
    }
}

void MatrixPrint(double mat[3][3]) {
    int i;
    printf("___________________________\r\n");
    for (i = 0; i < 3; i = i + 1) {
        printf(" |  %3.4f  | %3.4f  |  %3.4f |\r\n", (double) mat[i][0], (double) mat[i][1], (double) mat[i][2]); //this prints the matrix line by line
        printf("___________________________\r\n");

    }
}

void MatrixAdd(double mat1[3][3], double mat2[3][3], double result[3][3]) {
    double psum; // This is the template from which most of the other functions were created
    int j, k; // That's why all the commented out lines look the same
    psum = 0;
    for (j = 0; j < 3; j = j + 1) { //I like this because it's an easy way to hit every single entry in the matrix
        for (k = 0; k < 3; k = k + 1) {
            psum = (mat1[j][k] + mat2[j][k]);
            result[j][k] = psum;
            psum = 0;
        }
    }
}

int MatrixEquals(double mat1[3][3], double mat2[3][3]) {
    int j, k, ver;
    ver = 1;
    for (j = 0; j < 3; j = j + 1) { //if I wanted to make this work on a mxn matrix, all I'd need to do is change j to m and k to n
        for (k = 0; k < 3; k = k + 1) {
            if ((fabs(mat1[j][k] - mat2[j][k])) < FP_DELTA) {
                ver = ver * 1; //this is my easy way of making an AND function iterate over many loops, it needs 
            } else {
                ver = ver * 0; //one failure makes the whole thing a wash, hence a long AND
            }
        }
    }
    return (ver);
}

void MatrixScalarMultiply(double x, double mat[3][3], double result[3][3]) {
    double psum; // this is nearly identical to MatrixScalarAdd
    int j, k;
    psum = 0;
    for (j = 0; j < 3; j = j + 1) {
        for (k = 0; k < 3; k = k + 1) {
            psum = x * (mat[j][k]); // this simply multiplies every entry by x
            result[j][k] = psum;
            psum = 0;
        }
    }
}

void MatrixScalarAdd(double x, double mat[3][3], double result[3][3]) {
    double psum;
    int j, k;
    psum = 0;
    for (j = 0; j < 3; j = j + 1) {
        for (k = 0; k < 3; k = k + 1) {
            psum = x + (mat[j][k]);
            result[j][k] = psum;
            psum = 0;
        }
    }
}

double MatrixDeterminant(double mat[3][3]) {
    double det, psum1, psum2, psum3;
    det = 0;
    psum1 = (mat[0][0]*((mat[1][1] * mat[2][2])-(mat[1][2] * mat[2][1])));
    psum2 = (-1 * mat[0][1]*((mat[1][0] * mat[2][2])-(mat[2][0] * mat[1][2])));
    psum3 = (mat[0][2]*((mat[1][0] * mat[2][1])-(mat[2][0] * mat[1][1])));
    det = psum1 + psum2 + psum3;
    return (det);
}

void MatrixTranspose(double mat[3][3], double result[3][3]) {
    int j, k; // That's why all the commented out lines look the same
    for (j = 0; j < 3; j = j + 1) {
        for (k = 0; k < 3; k = k + 1) {
            result[j][k] = (mat[k][j]);
        }
    }
}

double MatrixTrace(double mat[3][3]) {
    int i;
    double result;
    result = 0;
    for (i = 0; i < 3; i = i + 1) {
        result = result + mat[i][i];
    }
    return (result);
}

void MatrixAdjugate(double mat[3][3], double result[3][3]) {
    double cofact[3][3];
    cofact[0][0] = ((mat[1][1] * mat[2][2])-(mat[1][2] * mat[2][1]));
    cofact[1][0] = ((mat[0][2] * mat[2][1])-(mat[0][1] * mat[2][2]));
    cofact[2][0] = ((mat[0][1] * mat[1][2])-(mat[0][2] * mat[1][1]));
    cofact[0][1] = ((mat[1][2] * mat[2][0])-(mat[1][0] * mat[2][2]));
    cofact[1][1] = ((mat[0][0] * mat[2][2])-(mat[0][2] * mat[2][0]));
    cofact[2][1] = ((mat[0][2] * mat[1][0])-(mat[0][0] * mat[1][2]));
    cofact[0][2] = ((mat[1][0] * mat[2][1])-(mat[1][1] * mat[2][0]));
    cofact[1][2] = ((mat[0][1] * mat[2][0])-(mat[0][0] * mat[2][1]));
    cofact[2][2] = ((mat[0][0] * mat[1][1])-(mat[0][1] * mat[1][0]));
    double det;
    det = MatrixDeterminant(mat);
    if (det == 0) {
        int j, k; // That's why all the commented out lines look the same
        for (j = 0; j < 3; j = j + 1) {
            for (k = 0; k < 3; k = k + 1) {
                result[j][k] = 0;
            }
        }
    } else {
        det = (1.0 / (det));
        MatrixScalarMultiply(det, cofact, result);
    }
}

void MatrixInverse(double mat[3][3], double result[3][3]) {
    double adj[3][3];
    MatrixAdjugate(mat, adj);
    MatrixTranspose(adj, result);
}

void MatrixCopy(double in[3][3], double out[3][3]){
    int j, k;
    for (j = 0; j < 3; j = j + 1) {
        for (k = 0; k < 3; k = k + 1) {
            out[j][k] = in[j][k];
        }
    }
}

void VectorCopy(double source[3], double dest[3]){
    int i;
    for(i=0;i<3;i++){
        dest[i]=source[i];
    }
}

void MatrixVectorProduct(double mat[3][3], double vect[3], double prod[3]){
    int i;
    for(i=0;i<3;i++){
        prod[i]=(vect[0]*mat[i][0])+(vect[1]*mat[i][1])+(vect[2]*mat[i][2]);
        //printf("(%f * %f) + (%f * %f) + (%f * %f) = %f",vect[i],mat[i][0],vect[i][1],mat[i][1],vect[i][2],mat[i][2],prod[i]);
    }
}

void VectorScale(double vect[3], double scale, double out[3]){
    int i;
    for(i=0;i<3;i++){
        out[i]=vect[i]*scale;
    }
}

void VectorAdd(double vect1[3], double vect2[3], double sum[3]){
    int i;
    for(i=0;i<3;i++){
        sum[i]=vect1[i]+vect2[i];
    }
}

void VectorPrint(double vect[3]){
    printf(" |  %3.4f  | %3.4f  |  %3.4f |\r\n", (double) vect[0], (double) vect[1], (double) vect[2]);
}


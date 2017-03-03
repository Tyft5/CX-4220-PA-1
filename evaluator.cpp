/*
 * CX 4220 / CSE 6220 Introduction to High Performance Computing
 *              Programming Assignment 1
 * 
 *  Serial polynomial evaluation algorithm function implementations goes here
 * 
 */

#include <math.h>

double poly_evaluator(const double x, const int n, const double* constants){
    //set the initial sum
    double sum = 0;
    //for each element
    for (int i = 0; i < n; i++) {
    	//add to the constant times x^i to the total sum
        sum += constants[i] * pow(x, i);
    }
    //return the sum
    return sum;
}

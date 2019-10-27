#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>
#include "ModelParameters.h"

bool CheckUnitarity(double Mh22, double MHA2, double theta, double vs, double d2)
{
    double lambda = Mh12*pow(cos(theta),2)/(2.0*vev2) + Mh22*pow(sin(theta),2)/(2.0*vev2);
    double delta2 = Mh12*sin(2.0*theta)/(vev*vs) - Mh22*sin(2.0*theta)/(vev*vs);
    gsl_matrix *ScatterMatrix = gsl_matrix_alloc(5,5); // 5x5
    gsl_vector *EigenValue = gsl_vector_alloc(5); // 5;
    gsl_matrix *EigenVec = gsl_matrix_alloc(5,5);
    gsl_eigen_symmv_workspace *work = gsl_eigen_symmv_alloc(5);
    double data[5][5] = {
        {4.0*lambda,             sqrt(2.0)*lambda, sqrt(2.0)*lambda, delta2/(2.0*sqrt(2.0)), delta2/(2.0*sqrt(2.0))},
        {sqrt(2.0)*lambda,       3.0*lambda,       lambda,           delta2/4.0,             delta2/4.0},
        {sqrt(2.0)*lambda,       lambda,           3.0*lambda,       delta2/4.0,             delta2/4.0},
        {delta2/(2.0*sqrt(2.0)), delta2/4.0,       delta2/4.0,       3.0*d2/4.0,             d2/4.0},
        {delta2/(2.0*sqrt(2.0)), delta2/4.0,       delta2/4.0,       d2/4.0,                 3.0*d2/4.0}
    };
    for (int i = 0; i < 5; ++i)
    {
        for (int j = 0; j < 5; ++j)
        {
            gsl_matrix_set(ScatterMatrix,i,j,data[i][j]);
        }
    }
    gsl_eigen_symmv(ScatterMatrix,EigenValue,EigenVec,work);
    gsl_eigen_symmv_free(work);
    gsl_eigen_symmv_sort(EigenValue, EigenVec, GSL_EIGEN_SORT_ABS_DESC);//sort from high to low according to the absolute value 
    double a0max = 1/16.0/Pi*abs(gsl_vector_get(EigenValue,0));
    gsl_matrix_free(ScatterMatrix);
    gsl_matrix_free(EigenVec);
    gsl_vector_free(EigenValue);
    return a0max<0.5;
}
bool CheckStability(double Mh22, double MHA2, double theta, double vs, double d2)
{
    double lambda = Mh12*pow(cos(theta),2)/(2.0*vev2) + Mh22*pow(sin(theta),2)/(2.0*vev2);
    double delta2 = Mh12*sin(2.0*theta)/(vev*vs) - Mh22*sin(2.0*theta)/(vev*vs);
    if (lambda<0)
    {
        return false;
    }
    if (d2<0)
    {
        return false;
    }
    if (delta2<0)
    {
        if (lambda*d2<delta2*delta2)
        {
            return false;
        }
    }
    return true;
}

/*----------------------------------------------------------------------------
 ADOL-C -- Automatic Differentiation by Overloading in C++
 File:     speelpenning.cpp
 Revision: $Id: speelpenning.cpp 299 2012-03-21 16:08:40Z kulshres $
 Contents: speelpennings example, described in the manual

 Copyright (c) Andrea Walther, Andreas Griewank, Andreas Kowarz, 
               Hristo Mitev, Sebastian Schlenkrich, Jean Utke, Olaf Vogel 
  
 This file is part of ADOL-C. This software is provided as open source.
 Any use, reproduction, or distribution of the software constitutes 
 recipient's acceptance of the terms of the accompanying license file.
 
---------------------------------------------------------------------------*/

/****************************************************************************/
/*                                                                 INCLUDES */
#include <adolc/adouble.h>            // use of active doubles
#include <adolc/drivers/drivers.h>    // use of "Easy to Use" drivers
#include <adolc/taping.h>             // use of taping

#include <iostream>
using namespace std;

#include <cstdlib>
#include <math.h>

#define DIM 2
#define DIMY 2

template <class type>
void my_function(const type x[DIM], type y[DIMY]) {

	// some nonlinear function of x
	y[0] = x[0] * x[1] * x[1];
	y[1] = x[1] * x[1] - x[0]*x[0];

}

double calc_sum_jac_diff(const double xp[DIM], double **jac_auto) {

	double err = 0.0;
	double jac_exact[DIMY][DIM];
	jac_exact[0][0] = xp[1]*xp[1];
	jac_exact[1][0] = -2*xp[0];
	jac_exact[0][1] = 2*xp[0]*xp[1];
	jac_exact[1][1] = 2*xp[1];
    for(int i = 0; i < DIMY; i++)
    	for(int j = 0; j < DIM; j++)
    		err += fabs(jac_exact[i][j]-jac_auto[i][j]);
    return err;
}

int main() {

    //size_t tape_stats[STAT_SIZE];

    double *xp = new double[DIM];
    double *yp = new double[DIMY];
    adouble *x = new adouble[DIM];
    adouble *y = new adouble[DIMY];

    for(int i = 0; i < DIM; i++)
        xp[i] = (i+1.0)/(2.0+i);           // some initialization

    trace_on(1);                         // tag = 1, keep = 0 by default
    for(int i = 0; i < DIM; i++) {
        x[i] <<= xp[i];                  // or  x <<= xp outside the loop
    }
    my_function(x,y);
    my_function(xp,yp);
    cout << "yp[0] = " << yp[0] << "\nyp[1] = " << yp[1] << endl;
    for(int i = 0; i < DIMY; i++) {
    	y[i] >>= yp[i];
    }
    trace_off(1);

    //tapestats(1,tape_stats);             // reading of tape statistics
    //cout<<"maxlive "<< tape_stats[NUM_MAX_LIVES]<<"\n";
    // ..... print other tape stats

    double **jac;
    jac = new double*[DIMY];
    for(int i = 0; i < DIMY; i++)
    	jac[i] = new double[DIM];

    jacobian(1,DIMY,DIM,xp,jac);

    double err_jac = calc_sum_jac_diff(xp,jac);

    cout << err_jac <<" error in jacobian \n";

    delete[] x,xp,y,yp;
    for(int i = 0; i < DIMY; i++)
    	delete[] jac[i];
    delete[] jac;
    return 0;
}

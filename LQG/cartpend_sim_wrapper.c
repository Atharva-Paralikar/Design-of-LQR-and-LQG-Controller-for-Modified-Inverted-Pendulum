
/*
 * Include Files
 *
 */
#if defined(MATLAB_MEX_FILE)
#include "tmwtypes.h"
#include "simstruc_types.h"
#else
#include "rtwtypes.h"
#endif



/* %%%-SFUNWIZ_wrapper_includes_Changes_BEGIN --- EDIT HERE TO _END */
#include <math.h>
/* %%%-SFUNWIZ_wrapper_includes_Changes_END --- EDIT HERE TO _BEGIN */
#define u_width 1
#define y_width 1

/*
 * Create external references here.  
 *
 */
/* %%%-SFUNWIZ_wrapper_externs_Changes_BEGIN --- EDIT HERE TO _END */
/* extern double func(double a); */
/* %%%-SFUNWIZ_wrapper_externs_Changes_END --- EDIT HERE TO _BEGIN */

/*
 * Output function
 *
 */
void cartpend_sim_Outputs_wrapper(const real_T *u0,
			real_T *y0,
			const real_T *xC)
{
/* %%%-SFUNWIZ_wrapper_Outputs_Changes_BEGIN --- EDIT HERE TO _END */
y0[0] = xC[0];
y0[1] = xC[1];
y0[2] = xC[2];
y0[3] = xC[3];
y0[4] = xC[4];
y0[5] = xC[5];
/* %%%-SFUNWIZ_wrapper_Outputs_Changes_END --- EDIT HERE TO _BEGIN */
}

/*
 * Derivatives function
 *
 */
void cartpend_sim_Derivatives_wrapper(const real_T *u0,
			real_T *y0,
			real_T *dx,
			real_T *xC)
{
/* %%%-SFUNWIZ_wrapper_Derivatives_Changes_BEGIN --- EDIT HERE TO _END */
double M=1000;
double m1=100;
double m2=100;
double l1=20;
double l2=10;
double g=9.81;

double F=u0[0];
dx[0] = xC[1];
dx[1]=(F-(g/2)*(m1*sin(2*xC[2])+m2*sin(2*xC[4]))-(m1*l1*(xC[3]*xC[3])*sin(xC[2]))-(m2*l2*(xC[5]*xC[5])*sin(xC[4])))/(M+m1*((sin(xC[2]))*(sin(xC[2])))+m2*((sin(xC[4]))*(sin(xC[4]))));
dx[2]= xC[3]; 
dx[3]= (dx[1]*cos(xC[2])-g*(sin(xC[2])))/l1; 
dx[4]= xC[5]; 
dx[5]= (dx[1]*cos(xC[4])-g*(sin(xC[4])))/l2;
/* %%%-SFUNWIZ_wrapper_Derivatives_Changes_END --- EDIT HERE TO _BEGIN */
}


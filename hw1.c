#include<stdlib.h>
#include<stdio.h>
#include <math.h>

double euler_method(double x_i, double t, double delta, double f(double t, double x));

double exp_ODE(double t, double x);

double recurse_method(double x_0, double t_0, double t_n, double delta, double f(double t, double x), 
    double method(double x_i, double t, double delta, double f(double t, double x)));

double midpoint_method(double x_i, double t, double delta, double f(double t, double x));

double rungekutta_method(double x_i, double t, double delta, double f(double t, double x));

double iterate_method(double x_0, double t_0, double t_n, double delta, double f(double t, double x),
    double method (double x_i, double t, double delta, double f(double t, double x)));

int main(double argc, char *argv[]){
  double x_0=1;
  double delta=1/pow(2,10);
  double t_0=0;
  double eulx_n;
  double midx_n;
  double rkmx_n;
  eulx_n=iterate_method(x_0,t_0,1, delta,exp_ODE, euler_method);
  midx_n=iterate_method(x_0,t_0,1, delta,exp_ODE, midpoint_method);
  rkmx_n=iterate_method(x_0,t_0,1,delta,exp_ODE,rungekutta_method);
  printf("%f\n",eulx_n);
  printf("%f\n",midx_n);
  printf("%f\n",rkmx_n);
  return 0;
}

double euler_method(double x_i, double t, double delta,double f(double t, double x)){
  double func_x=f(t,x_i);
  return x_i+delta*func_x;
}

double exp_ODE(double t, double x){
  return x;
}

double recurse_method(double x_0, double t_0, double t_n, double delta, double f(double t, double x), 
    double method(double x_i, double t, double delta, double f(double t, double x))){
 if (t_0<t_n){
   double x=method(x_0, t_0, delta, f);
   t_0+=delta;
   return(iterate_method(x,t_0,t_n,delta,f,method));
 }
 else
   return x_0;
}

double iterate_method(double x_0, double t_0, double t_n, double delta, double f(double t, double x),
    double method (double x_i, double t, double delta, double f(double t, double x))){
    double x_i=x_0;
    double t_i=t_0;
    while (t_i<t_n){
      x_i=method(x_i,t_i,delta,f);
      t_i+=delta;
    }
    return x_i;
    }
double midpoint_method(double x_i, double t, double delta, double f(double t, double x)){
  double k_temp=x_i+(delta/2)*f(t,x_i);
  return x_i+delta*f(t+delta/2,k_temp);
}

double rungekutta_method(double x_i, double t, double delta, double f(double t, double x)){
  double k_1=f(t,x_i);
  double k_2=f(t+delta/2,x_i+delta*k_1/2);
  double k_3=f(t+delta/2,x_i+delta*k_2/2);
  double k_4=f(t+delta,x_i+delta*k_3);
  return x_i+(delta*(k_1+2*k_2+2*k_3+k_4))/6;
}

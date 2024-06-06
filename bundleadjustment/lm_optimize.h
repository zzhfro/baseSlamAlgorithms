#ifndef LM_OPTIMIZE_H
#define LM_OPTIMIZE_H

#include <iostream>
#include <Eigen/Dense>
#include <random>
#include <assert.h>
#include "camera.h"

class lm_optimzie
{
//since  has many ways to optimize how to calculate jacobin this is just a framework 
//stop condition
// least square problem
struct option
{
double lm_mse_threshold=1e-16;
double lm_delta_threshold=1e-8;

//truth region radius
double lm_lamda = 1000;

//max iterations
int lm_max_iterations=300;
};
option lm_option;
inline void lm_set_param(double lm_mse_threshold,double lm_delta_threshold,double lm_lamda,double lm_max_iterations)
{
   this->lm_option.lm_mse_threshold=lm_mse_threshold;
   this->lm_option.lm_delta_threshold=lm_delta_threshold;
   this->lm_option.lm_lamda=lm_lamda;
   this->lm_option.lm_max_iterations=lm_max_iterations;
}

//variables to optimize
/**
 * ming(theta)=min||true_value-f(theta)||^2
 * 
 * x:the value to be optimized
 * delta_x: update x
 * f_x function output of x
 * 
 * 
*/
/*
Eigen::VectorXd x;
Eigen::VectorXd true_value;
Eigen::VectorXd delta_x;
Eigen::VectorXd f_x;
Eigen::VectorXd jacobin_x;
Eigen::MatrixXd JT_J;
*/


typedef void(*compute_f)(const Eigen::VectorXd &x,
                         Eigen::VectorXd &f_x);


typedef void(*compute_jacobin)(const Eigen::VectorXd &x,
                               Eigen::VectorXd &jaconbin);



typedef void(*solve_incremental_equation)(const Eigen::VectorXd &jacobin,const Eigen::VectorXd &error,
                                          Eigen::VectorXd &delta_x,
                                          double lm_lamda);

typedef void(*update_x)(Eigen::VectorXd x,const Eigen::VectorXd deltax);


void lm_process(compute_f f,compute_jacobin j_calculate,solve_incremental_equation incremental_calculate,update_x update,
                Eigen::VectorXd &error,Eigen::VectorXd &jacobin_x,Eigen::VectorXd &x,Eigen::VectorXd &delta_x )
{
  double lamda=1e-3;
  double error;
  int iterations=0;
  double error_new;
  Eigen::VectorXd fx_update;
  f(x,f_x);
  error=(true_value-f_x).squaredNorm();

  j_calculate(x,jacobin_x);

  incremental_calculate(jacobin_x,delta_x,x,f_x,lamda);
  while(lm_max_iterations>iterations)
  {  
    f(x+delta_x,fx_update);
     error_new=(true_value-fx_update).squaredNorm();
     if(error_new<error)
     {
      error=error_new;
      update(x,delta_x);
      if(delta_x.norm()<lm_delta_threshold)
      {
        break;
      }
      lamda=0.1*lamda;
      j_calculate(x,jacobin_x,f);
     }
     else
     {
      lamda=10*lamda;
     }
     incremental_calculate(jacobin_x,delta_x,x,f_x,lamda);
     iterations++;
  }
}










};



#endif 
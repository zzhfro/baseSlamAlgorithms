#ifndef BUNDLE_ADJUSTMENT_H
#define BUNDLE_ADJUSTMENT_H


#include <iostream>
#include <Eigen/Dense>
#include <random>
#include <assert.h>
#include "camera.h"
/**
 * a hand written BA
 * i use lm optmization method
 * achieve Jacobin calculate and lm optimization 
 * 
 * 
*/
/**
 *************************************************architecture description*******************************
 * 
 *struct bundle_adjustment_data
 *store the variable
 * 
 * 
 * 
 * camera_num  : int ,
 *               the number of camera
 * point3d_num : int ,
 *               the number of points
 * observations_num: int ,
 *                   the number of the observations point
 * 
 * point_w     : Eigen::MatrixXd    point3d_num,3 
 *               the point's 3d position in the world coordinate
 * 
 * observations_2d: Eigen::MatrixXd observations_num,2
 *                 pixel coordinates obtained by projecting a three-dimensional point onto the pixel coordinate system
 *                 the data organized in observations_2d in this way:
 *                 the point i observaed by camera j pij
 *                 p11,p21,.... pi2,pi+x,2,.....
 *            
 * 
 * 
 *point_index: Eigen::VectorXd observations,1
 *              store:index of the 3D point corresponding to the observation point
 * 
 *camera_index: Eigen::VectorXd observations,1
 *               store:index of the camera corresponding to the observation point
 *                     
 * 
 *Camera::ptr cameras: Vector<Camera> camera_num
 *                    
 *intermediate variable:
 * 
 * 
 *jacobian_x: Eigen::MatrixXd observations_num*2,camera_num*6+point_num*3
 * 
 *             this BA optimize the camera extrinsic parameters(2*6) and the 3d position of 3d point
 * 
 *x: Eigen::VectorXd camera_num*6+point_num*3,1
 *            
 *delta_x: Eigen::VectorXd camera_num*6+point_num*3,1
 * 
 *point2d_reprojected: Eigen::VectorXd observations*2,1
 *              the 3d point projected into the camera(reprojected)
 *              u1,v1,u2,v2,...,un,vn
 *               
 *point2d_raw: Eigen::VectorXd observations*2,1
 *                  the 3d point projected into the camera(observed)
 * 
 * 
 * ba_mode :indicate optimization type
 * 
 * 
 * 
 * 
 * BA,work flow:
 * 
 * adopt hand-writted lm method to perform non-linear optimize
 * 
 * step1: set the params of lm method(lm_set_param function)
 * 
 * step2: set the data of ba:   
 * 
 * step3: 
 * lm optimize begin:
 * it's work flow:
 *       call lm_process function
 *            the input params:   
 *                function ptr :
 *                    compute_f , input:  x: Eigen::VectorXd camera_num*6+point_num*3,1
 *                                output: projection2d_reprojected,Eigen::VectorXd 2*observations,1
 * 
 *                    compute_jacobin, input: input: x: Eigen::VectorXd camera_num*6+point_num*3,1
 *                                            output: jacobian_x,Eigen::matrixXd observations_num*2,camera_num*6+point_num*3
 * 
 *                    solve_incremental_equation: input:point2d_reprojected,point2d_raw,jacobin,lm_lamda,
 *                                                output:delta_x 
 *                    
 *                    update_x  : input:Eigen::VectorXd x, Eigen::VectorXd deltax
 *                    
 *                    
 * 
 *                                     
 *                                                       
 * 
*/
struct bundle_adjustment_data
{
  
  int camera_num;
  int point3d_num;
  int observations_num;

  Eigen::VectorXd camera_index;
  Eigen::VectorXd point_index;
  

  Eigen::MatrixXd point_w;

  Camera::camera_vector cameras;
  
  Eigen::MatrixXd jacobian_x;
  Eigen::VectorXd x;
  Eigen::VectorXd delta_x;

  Eigen::VectorXd point2d_raw;
  Eigen::VectorXd point2d_reprojected; 
  Eigen::VectorXd error;

};
struct option
{
double lm_mse_threshold=1e-16;
double lm_delta_threshold=1e-8;

//truth region radius
double lm_lamda= 1e-3;

//max iterations
int lm_max_iterations=300;
};



class bundle_adjustment
{
public:  
  enum ba_mode
  {
   BA_CAMERAS = 1,
   BA_POINTS  = 2,
   BA_CAMERAS_AND_POINTS = 3
  };
  bundle_adjustment_data data;
  option lm_option;



inline void set_data( int camera_num,int point3d_num,int observations_num,
                        Eigen::VectorXd &camera_index,Eigen::VectorXd &point_index,
                        Eigen::MatrixXd &point_w,Eigen::MatrixXd &observations_2d,
                        Camera::camera_vector &cameras,
                        Eigen::VectorXd &projection2d_raw
                      )
{
        this->data.camera_num=camera_num;
        this->data.point3d_num=point3d_num;
        this->data.observations_num=observations_num;

        this->data.camera_index=camera_index;
        this->data.point_index=point_index;

        this->data.point_w=point_w;
        this->data.cameras=cameras;

        this->data.point2d_raw.resize(2*this->data.observations_num);

        for(int i=0;i<2*observations_num;+i)
        {
          this->data.point2d_raw(2*i)=observations_2d(i,0);
          this->data.point2d_raw(2*i+1)=observations_2d(i,1);
          
        }

}

inline void lm_set_option(double lm_mse_threshold,double lm_delta_threshold,double lm_lamda,double lm_max_iterations)
{
        this->lm_option.lm_mse_threshold=lm_mse_threshold;
        this->lm_option.lm_delta_threshold=lm_delta_threshold;
        this->lm_option.lm_lamda=lm_lamda;
        this->lm_option.lm_max_iterations=lm_max_iterations;
}

inline void error_calculate();  

void jacobian_single_calculate(const Camera &c,const Eigen::Vector3d &point_w,
                                            Eigen::MatrixXd &camera_jacobian,Eigen::MatrixXd &point_jacobian);

void jacobian_calculate();    

void solve_incremental_equation(double lamda_inverse);

void update();

void ba_process();

                                       
  

};






#endif

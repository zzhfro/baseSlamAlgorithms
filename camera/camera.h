//created by ZH Zhao 2024.2

#ifndef CAMERA_H
#define CAMERA_H

#include<iostream>
#include <Eigen/Dense>
#include "sophus/se3.hpp"

class Camera{
      
public:
       typedef std::vector<Camera>  camera_vector;
       Camera()
       {
        c[0]=c[1]=0;
       }
       Eigen::Vector2d camera_projection(const Eigen::Vector3d &position_world) const
       {
         /**
          * get the position in camera coordinate
          * Rx'+t=x
          * x'=Rinverse(x-t)
          * */
         Eigen::Vector3d position_camera=camera_transform*position_world;
                                        
         /**
          * if consider the distortion
          * only consider Radial distortion
          * 
         */
         double x_hat=position_camera[0]/position_camera[2];
         double y_hat=position_camera[1]/position_camera[2]; 
         double r_square=x_hat*x_hat+y_hat*y_hat;
         double distortion_ratio = 1+ distortion_cofficient[0]*r_square+ distortion_cofficient[1]*r_square*r_square;
 
         /**
          * u=fx*x_dis+c[0]
          * v=fy*y_dis+c[1]
          * */ 
         Eigen::Vector2d u; 
         u[0]=f[0]*x_hat*distortion_ratio+c[0];
         u[1]=f[1]*y_hat*distortion_ratio+c[1];
         return u;
       }
       void intrinsic_set(double *f,double *c,double *distortion_cofficient)
       { 
          this->f[0]=f[0];
          this->f[1]=f[1];

          this->c[0]=c[0];
          this->c[1]=c[1];

          this->distortion_cofficient[0]=this->distortion_cofficient[0];
          this->distortion_cofficient[1]=this->distortion_cofficient[1];

       }
       // camera Intrinsics
       double f[2];
       double c[2];
       //Distortion coefficient
       double distortion_cofficient[2];
       //camera extrinsics
       //rotation translation is the projection of the camera

       Sophus::SE3d camera_transform;
       //the central coordinate
       
};
#endif
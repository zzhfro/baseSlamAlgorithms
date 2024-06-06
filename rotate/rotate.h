#ifndef ROTATE_H
#define ROTATE_H


#include <iostream>
#include <Eigen/Dense>
#include <assert.h>

//eigen has provided quaternion to axis_angle,rotaion to quaternion, quaternion rotate 

inline void axis_angle_rotate_point(const Eigen::Vector3d &rotation_axis,const Eigen::Vector3d &point,Eigen::Vector3d &point_rotate)
{
   double theta=rotation_axis.norm();
   Eigen::Matrix3d rotation_matrix;
   Eigen::Vector3d tmp;
   tmp=rotation_axis/theta;

   Eigen::MatrixXd I = Eigen::MatrixXd::Identity(3, 3);
   Eigen::MatrixXd n_delta;
   n_delta<<0,     -tmp[2], tmp[1],
            tmp[2], 0     ,-tmp[0],
            -tmp[2],tmp[0],0;
   
   double cos_theta=std::cos(theta);
   rotation_matrix=cos_theta*I+(1-cos_theta)*tmp*tmp.transpose()+std::sin(theta)*n_delta;
   point_rotate=rotation_matrix*point;

}

inline void quaternion_to_axis_angle(const Eigen::Quaterniond &q,Eigen::Vector3d &axis_angle)
{
  double theta=2*std::acos(q.w());
  double tmp=theta/std::sin(theta/2);
  axis_angle<<q.x()*tmp,q.y()*tmp,q.z()*tmp;
}



#endif

//created by ZH Zhao 2024.2
//triangulation , fundamental matrix estimation（using ransac）, get t and R from fundamental matrix
#ifndef EPIPOLAR_ALGORITHM_H
#define EPIPOLAR_ALGORITHM_H

#include <iostream>
#include <Eigen/Dense>
#include <random>
#include <assert.h>



void fundamental_matrix_calculate(Eigen::MatrixXd p1,Eigen::MatrixXd p2,Eigen::Matrix3d &fundamental_matrix);

double  sampson_distance(const Eigen::Vector2d &p1,const Eigen::Vector2d &p2,const Eigen::Matrix3d &f);

int  get_inliners_num(const Eigen::MatrixXd &p1,const Eigen::MatrixXd &p2,const Eigen::Matrix3d &f,double threshold,
                            Eigen::MatrixXd &inliner_p1,Eigen::MatrixXd &inliner_p2);

void random_sample(const Eigen::MatrixXd& p1,const Eigen::MatrixXd& p2, int num_rows,
                         Eigen::MatrixXd &p1_sample,Eigen::MatrixXd &p2_sample);

int fundamental_matrix_ransac(const Eigen::MatrixXd &p1,const Eigen::MatrixXd &p2,Eigen::Matrix3d &f); 

bool camera_pose_calculate(const Eigen::Matrix3d &f,const Eigen::Vector2d &p2d_c1,const Eigen::Vector2d &p2d_c2,
                           const Eigen::MatrixXd &k1,const Eigen::MatrixXd &k2, Eigen::Matrix3d &R, Eigen::Vector3d &t);                      

bool is_correct_pose(const Eigen::Vector2d &p2d_c1,const Eigen::Vector2d &p2d_c2,const Eigen::Matrix3d &R,const Eigen::Vector3d &t
                     ,const Eigen::Matrix3d &k1,const Eigen::Matrix3d &k2);
#endif
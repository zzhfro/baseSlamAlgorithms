#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Dense>
#include "epipolarAlgorithm.h"
using namespace Eigen;


int main()
{
 std::cout<<"step0"<<std::endl;
  Eigen::Vector2d p1={0.18012331426143646, -0.15658402442932129};
  Eigen::Vector2d p2={0.2082643061876297, -0.035404585301876068};

  Eigen::Matrix3d k1,k2;
  k1<<0.972222208,0,0,
      0,0.972222208,0,
      0,0,1;
  k2<<0.972222208,0,0,
      0,0.972222208,0,
      0,0,1;
      std::cout<<"step0"<<std::endl;
/**
 * test function
 * 
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
 * 
*/
Eigen::Matrix3d f;
std::cout<<"step0"<<std::endl;
f<<-0.005191866820221588,-0.015460923969578466,0.35260470328319654,
   0.022451443619913483,-0.00079225386526248181,-0.027885130552744289,
   -0.35188558059920161,0.032418724757766811,-0.005524537443406155;
Eigen::Matrix3d R;
Eigen::Vector3d t;
if(camera_pose_calculate(f,p1,p2,k1,k2,R,t))
{
  std::cout<<"R,t is "<<std::endl;
  std::cout<<R<<std::endl;
  std::cout<<t<<std::endl;
}




}   
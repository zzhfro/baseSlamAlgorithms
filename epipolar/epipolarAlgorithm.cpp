#include "epipolarAlgorithm.h"


/**formula and defination
 * three 3d position in 3coordinate p3d_w p3d_c1 p3d_c2
 * the 3d Point's position in world coordinate camera1 coordiante camera2 coordinate
 * 
 * 2d image(u,v): p2d_image1 p2d_image2
 * 
 * fundamental matrix:
 * p_c1=K[I,0]
 * p_c2=K[R,t]
 * E=t^R
 * F=K2-T*E*K1-1
 *     
 * triangulation:
 * int this file R,t for camera is projection matrix not the posture of camera
 * P_c=[R,t]*P_w
 * projection_matrix=K[R,t]
 *  
 * 
 * 
*/
void fundamental_matrix_calculate(Eigen::MatrixXd p1,Eigen::MatrixXd p2,Eigen::Matrix3d &fundamental_matrix)
{
  assert(p1.rows()==p2.rows());
  assert(p1.rows()>7);
  int point_num=p1.rows();
  Eigen::MatrixXd A(point_num, 9);
  //need to check
  for(int i=0;i<point_num;++i)
  {
    A(i, 0) = p1(i,0)*p2(i,0);
    A(i, 1) = p1(i,1)*p2(i,0);
    A(i, 2) = p2(i,0);
    A(i, 3) = p1(i,0)*p2(i,1);
    A(i, 4) = p1(i,1)*p2(i,1);
    A(i, 5) = p2(i,1);
    A(i, 6) = p1(i,0);
    A(i, 7) = p1(i,1);
    A(i, 8) = 1.0;
  }
  //SVD get the fundamental matrix
  //Af=0 A=UDVT f=vn
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeFullV);
  Eigen::MatrixXd V = svd.matrixV();
  Eigen::VectorXd f = V.col(V.cols() - 1);
  fundamental_matrix(0,0) = f[0]; fundamental_matrix(0,1) = f[1]; fundamental_matrix(0,2) = f[2];
  fundamental_matrix(1,0) = f[3]; fundamental_matrix(1,1) = f[4]; fundamental_matrix(1,2) = f[5];
  fundamental_matrix(2,0) = f[6]; fundamental_matrix(2,1) = f[7]; fundamental_matrix(2,2) = f[8];
  
  //apply singularity constraint (the third singularity should be 0)
  Eigen::JacobiSVD<Eigen::MatrixXd> svd_fundamental(fundamental_matrix, Eigen::ComputeFullU | Eigen::ComputeFullV);
  Eigen::MatrixXd U = svd_fundamental.matrixU();
  Eigen::MatrixXd V2;
  V2 = svd_fundamental.matrixV();
  Eigen::VectorXd singular_values = svd_fundamental.singularValues();
  singular_values(2) = 0;
  fundamental_matrix = U * singular_values.asDiagonal() * V2.transpose();
  


}

//ransac calculate fundamental_matrix

double  sampson_distance(const Eigen::Vector2d &p1,const Eigen::Vector2d &p2,const Eigen::Matrix3d &f)
{
   double distance=0;
   double p2_f_p1=p1[0]*(p2[0]*f(0,0)+p2[1]*f(1,0)+f(2,0))
                 +p1[1]*(p2[0]*f(0,1)+p2[1]*f(1,1)+f(2,1))
                        +p2[0]*f(0,2)+p2[1]*f(1,2)+f(2,2);
    p2_f_p1*=p2_f_p1;
    double tmp=0;
    tmp+=std::pow(f(0,0)*p1[0]+f(0,1)*p1[1]+f(0,2),2);
    tmp+=std::pow(f(1,0)*p1[0]+f(1,1)*p1[1]+f(1,2),2);
    tmp+=std::pow(f(0,0)*p2[0]+f(0,1)*p2[1]+f(0,2),2);
    tmp+=std::pow(f(1,0)*p2[0]+f(1,1)*p2[1]+f(1,2),2);
    
    return p2_f_p1/tmp;
   
}
int  get_inliners_num(const Eigen::MatrixXd &p1,const Eigen::MatrixXd &p2,const Eigen::Matrix3d &f,double threshold,
                            Eigen::MatrixXd &inliner_p1,Eigen::MatrixXd &inliner_p2)
{  
 int inliner_num=0;
 assert(p1.rows()==p2.rows());
 double error=0;
 for (int i=0;i<p1.rows();++i)
 {
  error=sampson_distance(p1.row(i),p2.row(i),f);
  if(error<threshold)
  {
   ++inliner_num;
   inliner_p1.row(i)=p1.row(i);
   inliner_p2.row(i)=p2.row(i);
  }
 }
 return inliner_num;

}

void random_sample(const Eigen::MatrixXd& p1,const Eigen::MatrixXd& p2, int num_rows,
                              Eigen::MatrixXd &p1_sample,Eigen::MatrixXd &p2_sample) 
{
    std::random_device rd;
    std::mt19937 gen(rd());

    Eigen::MatrixXd result(num_rows, p1.cols());

    std::vector<int> indices(p1.rows());
    std::iota(indices.begin(), indices.end(), 0); 

    std::shuffle(indices.begin(), indices.end(), gen); 
    std::cout<<"indices"<<std::endl;
    
    for (int i = 0; i < num_rows; ++i) {
     std::cout<<indices[i]<<std::endl;
        p1_sample.row(i) = p1.row(indices[i]); 
        p2_sample.row(i) = p2.row(indices[i]); 
    }
}
int fundamental_matrix_ransac(const Eigen::MatrixXd &p1,const Eigen::MatrixXd &p2,Eigen::Matrix3d &f)
{
 
  assert(p1.rows()==p2.rows());
  
  Eigen::MatrixXd eight_point_p1(8, 2);
  Eigen::MatrixXd eight_point_p2(8, 2);
  Eigen::MatrixXd inliner_p1_best(p1.rows(), 2);
  Eigen::MatrixXd inliner_p2_best(p1.rows(), 2);
  Eigen::MatrixXd inliner_p1_current(p1.rows(), 2);
  Eigen::MatrixXd inliner_p2_current(p1.rows(), 2);
  
  double p_inliner=0.5;
  double z=0.99;
  
  int sample_num=std::round(std::log(1-z)/std::log(1-std::pow(p_inliner,8)));
  int most_inliner_num=0;
  int current_inliner_num=0;
  double threshold=0.0015*0.0015;
  
 std::cout<<"step1"<<std::endl;
 
 for(int j=0;j<sample_num;++j)
 {
   
   random_sample(p1,p2,8,eight_point_p1,eight_point_p2);
   std::cout<<"random_sample"<<std::endl;
   std::cout<<eight_point_p1<<std::endl;
   fundamental_matrix_calculate(eight_point_p1,eight_point_p2,f);
   current_inliner_num=get_inliners_num(p1,p2,f,threshold,inliner_p1_current,inliner_p2_current);
   
    std::cout<<"sep0"<<std::endl;
    std::cout<<"current_inliner_num"<<std::endl;
    std::cout<<current_inliner_num<<std::endl;
   if(current_inliner_num>most_inliner_num)
   {
    most_inliner_num=current_inliner_num;
    inliner_p1_best=inliner_p1_current;
    inliner_p2_best=inliner_p2_current;
   }
   std::cout<<"sep1"<<std::endl;
 }
 std::cout<<"sep2"<<std::endl;
 std::cout<<inliner_p1_best.topRows(8)<<std::endl;
 std::cout<<"most_inliner_num"<<std::endl;
 std::cout<<most_inliner_num<<std::endl;
 fundamental_matrix_calculate(inliner_p1_best.topRows(most_inliner_num),inliner_p2_best.topRows(most_inliner_num),f);
 return most_inliner_num;
}

//分解E，F得到t，R，K

//triangulation (two point)
Eigen::Vector3d triangulation(const Eigen::Vector2d &p1,const Eigen::Vector2d &p2,
                              const Eigen::Matrix3d &R1_camera,const Eigen::Matrix3d &R2_camera,
                              const Eigen::Vector3d &t1_camera,const Eigen::Vector3d &t2_camera,
                              const Eigen::Matrix3d &k1,const Eigen::Matrix3d &k2)
                              
{


  Eigen::MatrixXd projection_1(3,4),projection_2(3,4);

  projection_1<<k1*R1_camera,k1*(t1_camera);
  projection_2<<k2*R2_camera,k2*(t2_camera);
  
  Eigen::MatrixXd A(4,4);
  A.row(0)=projection_1.row(2)*p1[0]-projection_1.row(0);
  A.row(1)=projection_1.row(2)*p1[1]-projection_1.row(1);
  A.row(2)=projection_2.row(2)*p2[0]-projection_2.row(0);
  A.row(3)=projection_2.row(2)*p2[1]-projection_2.row(1);
  
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeFullV);
  Eigen::MatrixXd V = svd.matrixV();
  Eigen::VectorXd P = V.col(V.cols() - 1);
  Eigen::Vector3d Pw;
  Pw<<P[0]/P[3],P[1]/P[3],P[2]/P[3];
  
  return Pw;
  
}    

bool camera_pose_calculate(const Eigen::Matrix3d &f,const Eigen::Vector2d &p2d_c1,const Eigen::Vector2d &p2d_c2,
                           const Eigen::MatrixXd &k1,const Eigen::MatrixXd &k2, Eigen::Matrix3d &R, Eigen::Vector3d &t)
{
  Eigen::Matrix3d k2_inverse=k2.inverse();
  
  Eigen::Matrix3d essential_matrix=k2_inverse.transpose()*f*k1.inverse();
  
  
  Eigen::JacobiSVD<Eigen::Matrix3d>svd(essential_matrix,Eigen::ComputeFullU | Eigen::ComputeFullV);
  
  Eigen::Matrix3d Rz;
  Rz<<0,-1,0,
      1,0,0,
      0,0,1;
  Eigen::Matrix3d Rz2;
  Rz2<<0,1,0
       -1,0,0,
       0,0,1;     
  Eigen::Matrix3d U=svd.matrixU();
  
  if(U.determinant()<0)
  {
    U.col(2)=-U.col(2);
  }
  Eigen::Matrix3d V=svd.matrixV();
  if(V.determinant()<0)
  {
    V.col(2)=-V.col(2);
  }
  
  Eigen::Vector3d t1=U.col(2);
  Eigen::Vector3d t2=-U.col(2);
  Eigen::Matrix3d R1=U*Rz*V.transpose();
  Eigen::Matrix3d R2=U*Rz2*V.transpose();
  int flags=0;
  
  if(is_correct_pose(p2d_c1,p2d_c2,R1,t1,k1,k2))
  {
   R=R1,t=t1;
   ++flags;
  }
  if(is_correct_pose(p2d_c1,p2d_c2,R1,t2,k1,k2))
  {
   R=R1,t=t2;
   ++flags;
  }
  if(is_correct_pose(p2d_c1,p2d_c2,R2,t1,k1,k2))
  {
   R=R2,t=t1;
   ++flags;
  }
  if(is_correct_pose(p2d_c1,p2d_c2,R2,t2,k1,k2))
  {
   R=R2,t=t2;
   ++flags;
  }
  return flags==1;
 
  
}

bool is_correct_pose(const Eigen::Vector2d &p2d_c1,const Eigen::Vector2d &p2d_c2,const Eigen::Matrix3d &R,const Eigen::Vector3d &t
                     ,const Eigen::Matrix3d &k1,const Eigen::Matrix3d &k2)
{
  //method : zc>0 for camera1 and camera2
  Eigen::MatrixXd R_c1 = Eigen::MatrixXd::Identity(3, 3);
  Eigen::Vector3d t_c1 = Eigen::Vector3d::Zero();
  Eigen::MatrixXd R_c2=R;
  Eigen::MatrixXd t_c2=t;
  Eigen::Vector3d P1;
  Eigen::Vector3d P2;
  
  P1=triangulation(p2d_c1,p2d_c2,R_c1,R_c2,t_c1,t_c2,k1,k2);//in fact it is P in camera1 coordinate
  
  P2=R_c2*P1+t_c2;
  return P1[2]>0.0f&&P2[2]>0.0f;


}

/**
 * 三角化要写ransac
 * PnP要写ransac
 * 
*/





#include "bundle_adjustment.h"
#include "lm_optimize.h"



inline void bundle_adjustment::error_calculate()
{
        Eigen::Vector2d tmp;
        for(int i=0;i<this->data.observations_num;++i)
        {
          tmp=this->data.cameras[this->data.camera_index(i)].camera_projection(this->data.point_w.row(this->data.point_index(i)));
          this->data.point2d_reprojected(2*i)=tmp(0);
          this->data.point2d_reprojected(2*i+1)=tmp(1);
          this->data.error(2*i)=this->data.point2d_raw(2*i)-tmp(0);
          this->data.error(2*i)=this->data.point2d_raw(2*i+1)-tmp(1);
        }
}



void bundle_adjustment::jacobian_single_calculate(const Camera &c,const Eigen::Vector3d &point_w,
                                            Eigen::MatrixXd &camera_jacobian,Eigen::MatrixXd &point_jacobian)
{
/**params 9 dims 
 * 0-2 quaternion
 * 3-5 translation
 * 7-9 focal length radial distortion  
 * 
 * currently we just optimize rotation and translations
 * point 3d location
 * 
 * 
*/          
            
          Eigen::VectorXd point_c=c.camera_transform*point_w;
          double f_x=c.f[0], f_y=c.f[1];

          double x=point_c[0],y=point_c[1],z=point_c[2];



          camera_jacobian<<-f_x/z,   0   , f_x*x/(z*z), f_x*x*y/(z*z)    , -(f_x+f_x*(x*x)/(z*z)) , f_x*y/z,
                          0,     -f_y/z , f_y*y/(z*z), f_y+f_y*y*y/(z*z), -f_y*x*y/(z*z)         , -f_y*x/z;


          Eigen::Matrix3d tmp;
          tmp<<-f_x/z,0       ,  f_x*x/(z*z),
              0     ,-f_y/z  ,  f_y*y/(z*z);
          point_jacobian=tmp*c.camera_transform.so3().matrix();     

} 


void bundle_adjustment::jacobian_calculate()
{
          Eigen::MatrixXd camera_jacobian(2*6);
          Eigen::MatrixXd point_jacobian(2*3);
          double camera_index,point_index;
          for(int i=0;i<this->data.observations_num;++i)
            { 
              camera_index=this->data.camera_index(i);
              point_index=this->data.point_index(i);
              jacobian_single_calculate(this->data.cameras[camera_index],this->data.point_w.row(point_index).transpose(),camera_jacobian,point_jacobian);
              this->data.jacobian_x.block(2*i,6*camera_index ,2,6) =camera_jacobian;
              this->data.jacobian_x.block(2*i,6*this->data.camera_num+3*point_index ,2,3) =point_jacobian;
            }
}

void bundle_adjustment::solve_incremental_equation(double lamda)   
{
      int point_num=this->data.point3d_num;
      int camera_num=this->data.camera_num;
      
      Eigen::MatrixXd jacobian_xx(3*point_num,3*point_num);
      Eigen::MatrixXd jacobian_cc(6*camera_num,6*camera_num);
      Eigen::MatrixXd jacobian_cx(6*camera_num,3*point_num);
      Eigen::MatrixXd jacobian_xc(3*point_num,6*camera_num);
      Eigen::MatrixXd jacobian_xx_inverse(3*point_num,3*point_num);

      Eigen::MatrixXd bc_coefficients(6*camera_num,6*camera_num);
      Eigen::MatrixXd bx_coefficients(3*point_num,3*point_num);

      Eigen::VectorXd bc(6*camera_num);
      Eigen::VectorXd bx(3*point_num);

      Eigen::VectorXd delta_camera(6*camera_num);
      Eigen::VectorXd delta_point(3*point_num);

      Eigen::VectorXd b;
      Eigen::MatrixXd H;
      Eigen::MatrixXd identity_matrix = Eigen::MatrixXd::Identity(6*camera_num+3*point_num, 6*camera_num+3*point_num);
      H=this->data.jacobian_x.transpose()*this->data.jacobian_x+lamda*identity_matrix;
      b=this->data.jacobian_x.transpose()*this->data.error;


      jacobian_xx=H.block(6*camera_num,6*camera_num,3*point_num,3*point_num);
      jacobian_cc=H.block(0,0,6*camera_num,6*camera_num);
      jacobian_cx=H.block(0,6*camera_num,6*camera_num,3*point_num);
      jacobian_xc=H.block(6*camera_num,0,3*point_num,6*camera_num);
      jacobian_xx_inverse=jacobian_xx.inverse();
      
      bc=b.segment(0,6*camera_num);
      bx=b.segment(6*camera_num,3*point_num);

      bc_coefficients=jacobian_cc-jacobian_cx*jacobian_xx_inverse*jacobian_xc;

      delta_camera= bc_coefficients.partialPivLu().solve(bc);
      delta_point= jacobian_xx_inverse*(bx-jacobian_xc*delta_camera);

      this->data.delta_x.segment(0,6*camera_num)=delta_camera;
      this->data.delta_x.segment(6*camera_num,3*point_num)=delta_point;


}

//update all the 3d point and camera's position   
void bundle_adjustment::update()
{
  
  Eigen::VectorXd camera_update(6);
  Eigen::VectorXd point_update(3);

  //update camera
  for(int i=0;i<this->data.camera_num;++i)
  {
    camera_update=this->data.delta_x.segment(6*i,6);
    this->data.cameras[i].camera_transform=Sophus::SE3d::exp(camera_update) * this->data.cameras[i].camera_transform;
  }
  
  
  //update point

  for(int i=0;i<this->data.point3d_num;++i)
  {
    this->data.point_w.row(i)=this->data.delta_x.segment(6*this->data.camera_num+i*3,3);
  }
}

void bundle_adjustment::ba_process()
{
  int iteration=0;
  double error=0;
  double error_new=0;
  
  Camera::camera_vector c_tmp_ptr=this->data.cameras;
  Eigen::MatrixXd point_w=this->data.point_w;

  error_calculate();
  error=this->data.error.squaredNorm();

  jacobian_calculate();
  solve_incremental_equation(this->lm_option.lm_lamda);

  while(this->lm_option.lm_max_iterations>iteration)
  { 
     update(); 
     error_calculate();
     error_new=this->data.error.squaredNorm();
     
     
     if(error_new<error)
     {
      c_tmp_ptr=this->data.cameras;
      point_w=this->data.point_w;
      error=error_new;
      if(this->data.delta_x.norm()<this->lm_option.lm_delta_threshold)
      {
        break;
      }
      this->lm_option.lm_lamda=0.1*this->lm_option.lm_lamda;
      jacobian_calculate();
     }
     
     else
     {
      this->data.cameras=c_tmp_ptr;
      this->data.point_w=point_w;
      this->lm_option.lm_lamda=10*this->lm_option.lm_lamda;
     }


     solve_incremental_equation(this->lm_option.lm_lamda);
     iteration++;
  }
}
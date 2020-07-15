#include "ros/ros.h"
#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/SVD"
#include "eigen3/Eigen/Core"

class KIN_LIB_IIWA {

   public:
      KIN_LIB_IIWA();
      ~KIN_LIB_IIWA();
      void compute_matrixA(double teta, int ith_link, Eigen::Matrix4d& A);
      void direct_kin(double q_curr[7], Eigen::Matrix4d& T);
      void compute_jacobian(double q_curr[7], Eigen::Matrix<double,6,7>& jacobian);
      void product_(Eigen::Vector3d f, Eigen::Vector4d c, Eigen::Vector4d d, Eigen::Vector3d& jp);
      void compute_w_dq(double teta[7],double w_d[7]);
   private:
      
     double a[7]; 
     double d[7];
     double alfa[7];
     // Eigen::Matrix4d A;
     // Eigen::Matrix4d T;
      //Eigen::Matrix<double,7,1> q_curr;
};



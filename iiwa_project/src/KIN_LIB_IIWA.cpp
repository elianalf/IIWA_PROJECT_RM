#include "KIN_LIB_IIWA.h"

using namespace std;

KIN_LIB_IIWA::KIN_LIB_IIWA(){
 
      alfa[0]= -1.57; //link 1. tra 1 e 0
      alfa[1]= 1.57; // link 2 
      alfa[2]= 1.57;
      alfa[3]= -1.57;
      alfa[4]= -1.57;
      alfa[5]= 1.57;
      alfa[6]= 0; //link 7
      
      
      d[0]=0.34; //link 1 - 0
      d[1]=0.0; //link 2
      d[2]=0.4; //link 3
      d[3]=0.0;// link 4
      d[4]=0.4; //link 5
      d[5]=0.0; //link6
      d[6]=0.126; //link 7, tra 7 e 6
      //d[6]=0.086;
      
      a[0]=0.0; //link 1
      a[1]= 0.0; 
      a[2]= 0.0;
      a[3]= 0.0;
      a[4]= 0.0;
      a[5]=0.0;
      a[6]=0.0;

 

}

KIN_LIB_IIWA::~KIN_LIB_IIWA(){}



void KIN_LIB_IIWA::compute_matrixA(double teta, int ith_link, Eigen::Matrix4d& A){
   
   A(0,0)=cos(teta);
   A(0,1)=-sin(teta)*cos(alfa[ith_link]);
   A(0,2)= sin(teta)*sin(alfa[ith_link]);
   A(0,3)=a[ith_link]*cos(teta); 
   
   A(1,0)=sin(teta); 
   A(1,1)=cos(teta)*cos(alfa[ith_link]);
   A(1,2) =-cos(teta)*sin(alfa[ith_link]); 
   A(1,3) =a[ith_link]*sin(teta);
    
   A(2,0)=0;
   A(2,1)= sin(alfa[ith_link]);
   A(2,2)= cos(alfa[ith_link]); 
   A(2,3)= d[ith_link];
    
   A(3,0)= 0; 
   A(3,1)= 0; 
   A(3,2)= 0; 
   A(3,3)= 1;

}


void KIN_LIB_IIWA::direct_kin(double q_curr[7], Eigen::Matrix4d& T){

   Eigen::Matrix4d A0; 
   Eigen::Matrix4d A1;
   Eigen::Matrix4d A2;
   Eigen::Matrix4d A3;
   Eigen::Matrix4d A4;
   Eigen::Matrix4d A5;
   Eigen::Matrix4d A6; 

   
   compute_matrixA(q_curr[0], 0,  A0);
   compute_matrixA(q_curr[1], 1,  A1);
   compute_matrixA(q_curr[2], 2,  A2);
   compute_matrixA(q_curr[3], 3,  A3);
   compute_matrixA(q_curr[4], 4,  A4);
   compute_matrixA(q_curr[5], 5,  A5);
   compute_matrixA(q_curr[6], 6,  A6);
     
   
   T = A0*A1*A2*A3*A4*A5*A6;
   

}


void KIN_LIB_IIWA::product_(Eigen::Vector3d f, Eigen::Vector4d c, Eigen::Vector4d d, Eigen::Vector3d& jp){
   Eigen::Vector3d b;
   for(int j=0; j<3; j++){
      b(j) =c(j)-d(j);
   }
   jp(0)=f(1)*b(2)-f(2)*b(1);
   jp(1)=f(2)*b(0)-f(0)*b(2);
   jp(2)=f(0)*b(1)-f(1)*b(0); 
}



void KIN_LIB_IIWA::compute_jacobian(double q_curr[7], Eigen::Matrix<double,6,7>& jacobian){
   //cout<<"in kin_lib"<<endl;
   
   Eigen::Matrix<double,3,1> z0;
   Eigen::Matrix<double,4,1> p0;
   Eigen::Matrix4d A0;
   Eigen::Matrix4d A1;
   Eigen::Matrix4d A2;
   Eigen::Matrix4d A3;
   Eigen::Matrix4d A4;
   Eigen::Matrix4d A5;
   Eigen::Matrix4d A6;
   z0<< 0, 0, 1;
   p0<< 0,0, 0, 1;
   // cout<<"-1"<<endl;
    
   compute_matrixA(q_curr[0], 0,  A0);
   compute_matrixA(q_curr[1], 1,  A1);
   compute_matrixA(q_curr[2], 2,  A2);
   compute_matrixA(q_curr[3], 3,  A3);
   compute_matrixA(q_curr[4], 4,  A4);
   compute_matrixA(q_curr[5], 5,  A5);
   compute_matrixA(q_curr[6], 6,  A6);
   //cout<<"-k"<<endl;
    
   Eigen::Matrix<double,3,1> z1;
   Eigen::Vector3d z2;
   Eigen::Vector3d z3;
   Eigen::Vector3d z4;
   Eigen::Vector3d z5;
   Eigen::Vector3d z6;
   Eigen::Vector3d z7;
 // cout<<"begin block"<<endl;
  
  
 
  z1=  A0.block<3,3>(0,0)*z0;
  z2=A0.block<3,3>(0,0)*A1.block<3,3>(0,0)*z0;
  z3=A0.block<3,3>(0,0)*A1.block<3,3>(0,0)*A2.block<3,3>(0,0)*z0;
  z4=A0.block<3,3>(0,0)*A1.block<3,3>(0,0)*A2.block<3,3>(0,0)*A3.block<3,3>(0,0)*z0;
  z5=A0.block<3,3>(0,0)*A1.block<3,3>(0,0)*A2.block<3,3>(0,0)*A3.block<3,3>(0,0)*A4.block<3,3>(0,0)*z0;
  z6=A0.block<3,3>(0,0)*A1.block<3,3>(0,0)*A2.block<3,3>(0,0)*A3.block<3,3>(0,0)*A4.block<3,3>(0,0)*A5.block<3,3>(0,0)*z0;
  z7=A0.block<3,3>(0,0)*A1.block<3,3>(0,0)*A2.block<3,3>(0,0)*A3.block<3,3>(0,0)*A4.block<3,3>(0,0)*A5.block<3,3>(0,0)*A6.block<3,3>(0,0)*z0;
  //cout<<"0"<<endl;
  
  Eigen::Vector4d p1=A0*p0;
  Eigen::Vector4d p2=A0*A1*p0;
  Eigen::Vector4d p3=A0*A1*A2*p0;
  Eigen::Vector4d p4=A0*A1*A2*A3*p0;
  Eigen::Vector4d p5=A0*A1*A2*A3*A4*p0;
  Eigen::Vector4d p6=A0*A1*A2*A3*A4*A5*p0;
  Eigen::Vector4d p7=A0*A1*A2*A3*A4*A5*A6*p0;
   //cout<<"1"<<endl;
   
   Eigen::Vector3d jp1;
   Eigen::Vector3d jp2;
   Eigen::Vector3d jp3;
   Eigen::Vector3d jp4;
   Eigen::Vector3d jp5;
   Eigen::Vector3d jp6;
   Eigen::Vector3d jp7;
  
   product_(z0, p7, p0, jp1);
   product_(z1, p7,p1, jp2);
   product_(z2, p7,p2, jp3);
   product_(z3, p7,p3, jp4);
   product_(z4, p7,p4, jp5);
   product_(z5, p7,p5, jp6);
   product_(z6, p7,p6, jp7);
   //cout<<"2"<<p4<<endl;
   Eigen::Matrix<double,6,1> J1;
   Eigen::Matrix<double,6,1> J2;
   Eigen::Matrix<double,6,1> J3;
   Eigen::Matrix<double,6,1> J4;
   Eigen::Matrix<double,6,1> J5;
   Eigen::Matrix<double,6,1> J6;
   Eigen::Matrix<double,6,1> J7;
   
   J1.topRows<3>() = jp1;   
   J2.topRows<3>() = jp2;
   J3.topRows<3>() = jp3;
   J4.topRows<3>() = jp4;
   J5.topRows<3>() = jp5;
   J6.topRows<3>() = jp6;
   J7.topRows<3>() = jp7;
   
   J1.bottomRows<3>() = z0;
   J2.bottomRows<3>() = z1;
   J3.bottomRows<3>() = z2;
   J4.bottomRows<3>() = z3;
   J5.bottomRows<3>() = z4;
   J6.bottomRows<3>() = z5;
   J7.bottomRows<3>() = z6;
   //cout<<"3"<<jp1<<endl;
   jacobian.block<6,1>(0,0) = J1;
   jacobian.block<6,1>(0,1) = J2;
   jacobian.block<6,1>(0,2) = J3;
   jacobian.block<6,1>(0,3) = J4;
   jacobian.block<6,1>(0,4) = J5;
   jacobian.block<6,1>(0,5) = J6;
   jacobian.block<6,1>(0,6) = J7;
}



void KIN_LIB_IIWA::compute_w_dq(double teta[7], double w_d[7]){

w_d[0] = -(6000*cos(teta[0])*sin(teta[1]) - 6000*sin(teta[0])*sin(teta[1]) - 6000*cos(teta[3])*sin(teta[0])*sin(teta[1]) + 6000*cos(teta[0])*sin(teta[2])*sin(teta[3]) + 6000*sin(teta[0])*sin(teta[2])*sin(teta[3]) + 6000*cos(teta[0])*cos(teta[3])*sin(teta[1]) - 6000*cos(teta[0])*cos(teta[1])*cos(teta[2])*sin(teta[3]) + 12900*cos(teta[0])*cos(teta[3])*cos(teta[5])*sin(teta[1]) + 6000*cos(teta[1])*cos(teta[2])*sin(teta[0])*sin(teta[3]) - 12900*cos(teta[3])*cos(teta[5])*sin(teta[0])*sin(teta[1]) + 12900*cos(teta[0])*cos(teta[5])*sin(teta[2])*sin(teta[3]) - 12900*cos(teta[0])*cos(teta[2])*sin(teta[4])*sin(teta[5]) + 12900*cos(teta[5])*sin(teta[0])*sin(teta[2])*sin(teta[3]) - 12900*cos(teta[2])*sin(teta[0])*sin(teta[4])*sin(teta[5]) - 12900*cos(teta[0])*cos(teta[1])*cos(teta[2])*cos(teta[5])*sin(teta[3]) + 12900*cos(teta[1])*cos(teta[2])*cos(teta[5])*sin(teta[0])*sin(teta[3]) - 12900*cos(teta[0])*cos(teta[3])*cos(teta[4])*sin(teta[2])*sin(teta[5]) - 12900*cos(teta[0])*cos(teta[1])*sin(teta[2])*sin(teta[4])*sin(teta[5]) + 12900*cos(teta[0])*cos(teta[4])*sin(teta[1])*sin(teta[3])*sin(teta[5]) - 12900*cos(teta[3])*cos(teta[4])*sin(teta[0])*sin(teta[2])*sin(teta[5]) + 12900*cos(teta[1])*sin(teta[0])*sin(teta[2])*sin(teta[4])*sin(teta[5]) - 12900*cos(teta[4])*sin(teta[0])*sin(teta[1])*sin(teta[3])*sin(teta[5]) + 12900*cos(teta[0])*cos(teta[1])*cos(teta[2])*cos(teta[3])*cos(teta[4])*sin(teta[5]) - 12900*cos(teta[1])*cos(teta[2])*cos(teta[3])*cos(teta[4])*sin(teta[0])*sin(teta[5]))/(200* sqrt( pow(40*sin(teta[0])*sin(teta[1]) - 86*cos(teta[5])*(cos(teta[0])*sin(teta[2])*sin(teta[3]) - cos(teta[3])*sin(teta[0])*sin(teta[1]) + cos(teta[1])*cos(teta[2])*sin(teta[0])*sin(teta[3])) - 40*sin(teta[3])*(cos(teta[0])*sin(teta[2]) + cos(teta[1])*cos(teta[2])*sin(teta[0])) + 86*sin(teta[5])*(cos(teta[4])*(sin(teta[0])*sin(teta[1])*sin(teta[3]) + cos(teta[0])*cos(teta[3])*sin(teta[2]) + cos(teta[1])*cos(teta[2])*cos(teta[3])*sin(teta[0])) + sin(teta[4])*(cos(teta[0])*cos(teta[2]) - cos(teta[1])*sin(teta[0])*sin(teta[2]))) + 40*cos(teta[3])*sin(teta[0])*sin(teta[1]) - 75, 2) + 4*  pow(20*cos(teta[1]) + 20*cos(teta[1])*cos(teta[3]) + 20*cos(teta[2])*sin(teta[1])*sin(teta[3]) + 43*cos(teta[1])*cos(teta[3])*cos(teta[5]) + 43*cos(teta[2])*cos(teta[5])*sin(teta[1])*sin(teta[3]) + 43*cos(teta[1])*cos(teta[4])*sin(teta[3])*sin(teta[5]) + 43*sin(teta[1])*sin(teta[2])*sin(teta[4])*sin(teta[5]) - 43*cos(teta[2])*cos(teta[3])*cos(teta[4])*sin(teta[1])*sin(teta[5]) + 7, 2) + pow(40*cos(teta[0])*sin(teta[1]) + 86*cos(teta[5])*(sin(teta[0])*sin(teta[2])*sin(teta[3]) + cos(teta[0])*cos(teta[3])*sin(teta[1]) - cos(teta[0])*cos(teta[1])*cos(teta[2])*sin(teta[3])) + 40*sin(teta[3])*(sin(teta[0])*sin(teta[2]) - cos(teta[0])*cos(teta[1])*cos(teta[2])) + 86*sin(teta[5])*(cos(teta[4])*(cos(teta[0])*sin(teta[1])*sin(teta[3]) - cos(teta[3])*sin(teta[0])*sin(teta[2]) + cos(teta[0])*cos(teta[1])*cos(teta[2])*cos(teta[3])) - sin(teta[4])*(cos(teta[2])*sin(teta[0]) + cos(teta[0])*cos(teta[1])*sin(teta[2]))) + 40*cos(teta[0])*cos(teta[3])*sin(teta[1]) - 75, 2)));


w_d[1]=(4*cos(teta[0])*(40*cos(teta[0])*sin(teta[1]) + 86*cos(teta[5])*(sin(teta[0])*sin(teta[2])*sin(teta[3]) + cos(teta[0])*cos(teta[3])*sin(teta[1]) - cos(teta[0])*cos(teta[1])*cos(teta[2])*sin(teta[3])) + 40*sin(teta[3])*(sin(teta[0])*sin(teta[2]) - cos(teta[0])*cos(teta[1])*cos(teta[2])) + 86*sin(teta[5])*(cos(teta[4])*(cos(teta[0])*sin(teta[1])*sin(teta[3]) - cos(teta[3])*sin(teta[0])*sin(teta[2]) + cos(teta[0])*cos(teta[1])*cos(teta[2])*cos(teta[3])) - sin(teta[4])*(cos(teta[2])*sin(teta[0]) + cos(teta[0])*cos(teta[1])*sin(teta[2]))) + 40*cos(teta[0])*cos(teta[3])*sin(teta[1]) - 75)*(20*cos(teta[1]) + 20*cos(teta[1])*cos(teta[3]) + 20*cos(teta[2])*sin(teta[1])*sin(teta[3]) + 43*cos(teta[1])*cos(teta[3])*cos(teta[5]) + 43*cos(teta[2])*cos(teta[5])*sin(teta[1])*sin(teta[3]) + 43*cos(teta[1])*cos(teta[4])*sin(teta[3])*sin(teta[5]) + 43*sin(teta[1])*sin(teta[2])*sin(teta[4])*sin(teta[5]) - 43*cos(teta[2])*cos(teta[3])*cos(teta[4])*sin(teta[1])*sin(teta[5])) - 8*(20*sin(teta[1]) + 20*cos(teta[3])*sin(teta[1]) - 20*cos(teta[1])*cos(teta[2])*sin(teta[3]) + 43*cos(teta[3])*cos(teta[5])*sin(teta[1]) - 43*cos(teta[1])*cos(teta[2])*cos(teta[5])*sin(teta[3]) - 43*cos(teta[1])*sin(teta[2])*sin(teta[4])*sin(teta[5]) + 43*cos(teta[4])*sin(teta[1])*sin(teta[3])*sin(teta[5]) + 43*cos(teta[1])*cos(teta[2])*cos(teta[3])*cos(teta[4])*sin(teta[5]))*(20*cos(teta[1]) + 20*cos(teta[1])*cos(teta[3]) + 20*cos(teta[2])*sin(teta[1])*sin(teta[3]) + 43*cos(teta[1])*cos(teta[3])*cos(teta[5]) + 43*cos(teta[2])*cos(teta[5])*sin(teta[1])*sin(teta[3]) + 43*cos(teta[1])*cos(teta[4])*sin(teta[3])*sin(teta[5]) + 43*sin(teta[1])*sin(teta[2])*sin(teta[4])*sin(teta[5]) - 43*cos(teta[2])*cos(teta[3])*cos(teta[4])*sin(teta[1])*sin(teta[5]) + 7) + 4*sin(teta[0])*(40*sin(teta[0])*sin(teta[1]) - 86*cos(teta[5])*(cos(teta[0])*sin(teta[2])*sin(teta[3]) - cos(teta[3])*sin(teta[0])*sin(teta[1]) + cos(teta[1])*cos(teta[2])*sin(teta[0])*sin(teta[3])) - 40*sin(teta[3])*(cos(teta[0])*sin(teta[2]) + cos(teta[1])*cos(teta[2])*sin(teta[0])) + 86*sin(teta[5])*(cos(teta[4])*(sin(teta[0])*sin(teta[1])*sin(teta[3]) + cos(teta[0])*cos(teta[3])*sin(teta[2]) + cos(teta[1])*cos(teta[2])*cos(teta[3])*sin(teta[0])) + sin(teta[4])*(cos(teta[0])*cos(teta[2]) - cos(teta[1])*sin(teta[0])*sin(teta[2]))) + 40*cos(teta[3])*sin(teta[0])*sin(teta[1]) - 75)*(20*cos(teta[1]) + 20*cos(teta[1])*cos(teta[3]) + 20*cos(teta[2])*sin(teta[1])*sin(teta[3]) + 43*cos(teta[1])*cos(teta[3])*cos(teta[5]) + 43*cos(teta[2])*cos(teta[5])*sin(teta[1])*sin(teta[3]) + 43*cos(teta[1])*cos(teta[4])*sin(teta[3])*sin(teta[5]) + 43*sin(teta[1])*sin(teta[2])*sin(teta[4])*sin(teta[5]) - 43*cos(teta[2])*cos(teta[3])*cos(teta[4])*sin(teta[1])*sin(teta[5])))/(200* sqrt( pow(40*sin(teta[0])*sin(teta[1]) - 86*cos(teta[5])*(cos(teta[0])*sin(teta[2])*sin(teta[3]) - cos(teta[3])*sin(teta[0])*sin(teta[1]) + cos(teta[1])*cos(teta[2])*sin(teta[0])*sin(teta[3])) - 40*sin(teta[3])*(cos(teta[0])*sin(teta[2]) + cos(teta[1])*cos(teta[2])*sin(teta[0])) + 86*sin(teta[5])*(cos(teta[4])*(sin(teta[0])*sin(teta[1])*sin(teta[3]) + cos(teta[0])*cos(teta[3])*sin(teta[2]) + cos(teta[1])*cos(teta[2])*cos(teta[3])*sin(teta[0])) + sin(teta[4])*(cos(teta[0])*cos(teta[2]) - cos(teta[1])*sin(teta[0])*sin(teta[2]))) + 40*cos(teta[3])*sin(teta[0])*sin(teta[1]) - 75, 2) + 4* pow(20*cos(teta[1]) + 20*cos(teta[1])*cos(teta[3]) + 20*cos(teta[2])*sin(teta[1])*sin(teta[3]) + 43*cos(teta[1])*cos(teta[3])*cos(teta[5]) + 43*cos(teta[2])*cos(teta[5])*sin(teta[1])*sin(teta[3]) + 43*cos(teta[1])*cos(teta[4])*sin(teta[3])*sin(teta[5]) + 43*sin(teta[1])*sin(teta[2])*sin(teta[4])*sin(teta[5]) - 43*cos(teta[2])*cos(teta[3])*cos(teta[4])*sin(teta[1])*sin(teta[5]) + 7 , 2) + pow(40*cos(teta[0])*sin(teta[1]) + 86*cos(teta[5])*(sin(teta[0])*sin(teta[2])*sin(teta[3]) + cos(teta[0])*cos(teta[3])*sin(teta[1]) - cos(teta[0])*cos(teta[1])*cos(teta[2])*sin(teta[3])) + 40*sin(teta[3])*(sin(teta[0])*sin(teta[2]) - cos(teta[0])*cos(teta[1])*cos(teta[2])) + 86*sin(teta[5])*(cos(teta[4])*(cos(teta[0])*sin(teta[1])*sin(teta[3]) - cos(teta[3])*sin(teta[0])*sin(teta[2]) + cos(teta[0])*cos(teta[1])*cos(teta[2])*cos(teta[3])) - sin(teta[4])*(cos(teta[2])*sin(teta[0]) + cos(teta[0])*cos(teta[1])*sin(teta[2]))) + 40*cos(teta[0])*cos(teta[3])*sin(teta[1]) - 75,2)  ));




w_d[2]=-(2*(40*sin(teta[3])*(cos(teta[0])*cos(teta[2]) - cos(teta[1])*sin(teta[0])*sin(teta[2])) + 86*sin(teta[5])*(sin(teta[4])*(cos(teta[0])*sin(teta[2]) + cos(teta[1])*cos(teta[2])*sin(teta[0])) - cos(teta[3])*cos(teta[4])*(cos(teta[0])*cos(teta[2]) - cos(teta[1])*sin(teta[0])*sin(teta[2]))) + 86*cos(teta[5])*sin(teta[3])*(cos(teta[0])*cos(teta[2]) - cos(teta[1])*sin(teta[0])*sin(teta[2])))*(40*sin(teta[0])*sin(teta[1]) - 86*cos(teta[5])*(cos(teta[0])*sin(teta[2])*sin(teta[3]) - cos(teta[3])*sin(teta[0])*sin(teta[1]) + cos(teta[1])*cos(teta[2])*sin(teta[0])*sin(teta[3])) - 40*sin(teta[3])*(cos(teta[0])*sin(teta[2]) + cos(teta[1])*cos(teta[2])*sin(teta[0])) + 86*sin(teta[5])*(cos(teta[4])*(sin(teta[0])*sin(teta[1])*sin(teta[3]) + cos(teta[0])*cos(teta[3])*sin(teta[2]) + cos(teta[1])*cos(teta[2])*cos(teta[3])*sin(teta[0])) + sin(teta[4])*(cos(teta[0])*cos(teta[2]) - cos(teta[1])*sin(teta[0])*sin(teta[2]))) + 40*cos(teta[3])*sin(teta[0])*sin(teta[1]) - 75) - 2*(40*sin(teta[3])*(cos(teta[2])*sin(teta[0]) + cos(teta[0])*cos(teta[1])*sin(teta[2])) - 86*sin(teta[5])*(cos(teta[0])*cos(teta[1])*cos(teta[2])*sin(teta[4]) - sin(teta[0])*sin(teta[2])*sin(teta[4]) + cos(teta[2])*cos(teta[3])*cos(teta[4])*sin(teta[0]) + cos(teta[0])*cos(teta[1])*cos(teta[3])*cos(teta[4])*sin(teta[2])) + 86*cos(teta[5])*sin(teta[3])*(cos(teta[2])*sin(teta[0]) + cos(teta[0])*cos(teta[1])*sin(teta[2])))*(40*cos(teta[0])*sin(teta[1]) + 86*cos(teta[5])*(sin(teta[0])*sin(teta[2])*sin(teta[3]) + cos(teta[0])*cos(teta[3])*sin(teta[1]) - cos(teta[0])*cos(teta[1])*cos(teta[2])*sin(teta[3])) + 40*sin(teta[3])*(sin(teta[0])*sin(teta[2]) - cos(teta[0])*cos(teta[1])*cos(teta[2])) + 86*sin(teta[5])*(cos(teta[4])*(cos(teta[0])*sin(teta[1])*sin(teta[3]) - cos(teta[3])*sin(teta[0])*sin(teta[2]) + cos(teta[0])*cos(teta[1])*cos(teta[2])*cos(teta[3])) - sin(teta[4])*(cos(teta[2])*sin(teta[0]) + cos(teta[0])*cos(teta[1])*sin(teta[2]))) + 40*cos(teta[0])*cos(teta[3])*sin(teta[1]) - 75) + 8*sin(teta[1])*(20*sin(teta[2])*sin(teta[3]) + 43*cos(teta[5])*sin(teta[2])*sin(teta[3]) - 43*cos(teta[2])*sin(teta[4])*sin(teta[5]) - 43*cos(teta[3])*cos(teta[4])*sin(teta[2])*sin(teta[5]))*(20*cos(teta[1]) + 20*cos(teta[1])*cos(teta[3]) + 20*cos(teta[2])*sin(teta[1])*sin(teta[3]) + 43*cos(teta[1])*cos(teta[3])*cos(teta[5]) + 43*cos(teta[2])*cos(teta[5])*sin(teta[1])*sin(teta[3]) + 43*cos(teta[1])*cos(teta[4])*sin(teta[3])*sin(teta[5]) + 43*sin(teta[1])*sin(teta[2])*sin(teta[4])*sin(teta[5]) - 43*cos(teta[2])*cos(teta[3])*cos(teta[4])*sin(teta[1])*sin(teta[5]) + 7))/(200* sqrt( pow(40*sin(teta[0])*sin(teta[1]) - 86*cos(teta[5])*(cos(teta[0])*sin(teta[2])*sin(teta[3]) - cos(teta[3])*sin(teta[0])*sin(teta[1]) + cos(teta[1])*cos(teta[2])*sin(teta[0])*sin(teta[3])) - 40*sin(teta[3])*(cos(teta[0])*sin(teta[2]) + cos(teta[1])*cos(teta[2])*sin(teta[0])) + 86*sin(teta[5])*(cos(teta[4])*(sin(teta[0])*sin(teta[1])*sin(teta[3]) + cos(teta[0])*cos(teta[3])*sin(teta[2]) + cos(teta[1])*cos(teta[2])*cos(teta[3])*sin(teta[0])) + sin(teta[4])*(cos(teta[0])*cos(teta[2]) - cos(teta[1])*sin(teta[0])*sin(teta[2]))) + 40*cos(teta[3])*sin(teta[0])*sin(teta[1]) - 75, 2) + 4* pow(20*cos(teta[1]) + 20*cos(teta[1])*cos(teta[3]) + 20*cos(teta[2])*sin(teta[1])*sin(teta[3]) + 43*cos(teta[1])*cos(teta[3])*cos(teta[5]) + 43*cos(teta[2])*cos(teta[5])*sin(teta[1])*sin(teta[3]) + 43*cos(teta[1])*cos(teta[4])*sin(teta[3])*sin(teta[5]) + 43*sin(teta[1])*sin(teta[2])*sin(teta[4])*sin(teta[5]) - 43*cos(teta[2])*cos(teta[3])*cos(teta[4])*sin(teta[1])*sin(teta[5]) + 7, 2) + pow(40*cos(teta[0])*sin(teta[1]) + 86*cos(teta[5])*(sin(teta[0])*sin(teta[2])*sin(teta[3]) + cos(teta[0])*cos(teta[3])*sin(teta[1]) - cos(teta[0])*cos(teta[1])*cos(teta[2])*sin(teta[3])) + 40*sin(teta[3])*(sin(teta[0])*sin(teta[2]) - cos(teta[0])*cos(teta[1])*cos(teta[2])) + 86*sin(teta[5])*(cos(teta[4])*(cos(teta[0])*sin(teta[1])*sin(teta[3]) - cos(teta[3])*sin(teta[0])*sin(teta[2]) + cos(teta[0])*cos(teta[1])*cos(teta[2])*cos(teta[3])) - sin(teta[4])*(cos(teta[2])*sin(teta[0]) + cos(teta[0])*cos(teta[1])*sin(teta[2]))) + 40*cos(teta[0])*cos(teta[3])*sin(teta[1]) - 75, 2 )) );


w_d[3]= -(2*(86*cos(teta[5])*(sin(teta[0])*sin(teta[1])*sin(teta[3]) + cos(teta[0])*cos(teta[3])*sin(teta[2]) + cos(teta[1])*cos(teta[2])*cos(teta[3])*sin(teta[0])) + 40*cos(teta[3])*(cos(teta[0])*sin(teta[2]) + cos(teta[1])*cos(teta[2])*sin(teta[0])) + 40*sin(teta[0])*sin(teta[1])*sin(teta[3]) + 86*cos(teta[4])*sin(teta[5])*(cos(teta[0])*sin(teta[2])*sin(teta[3]) - cos(teta[3])*sin(teta[0])*sin(teta[1]) + cos(teta[1])*cos(teta[2])*sin(teta[0])*sin(teta[3])))*(40*sin(teta[0])*sin(teta[1]) - 86*cos(teta[5])*(cos(teta[0])*sin(teta[2])*sin(teta[3]) - cos(teta[3])*sin(teta[0])*sin(teta[1]) + cos(teta[1])*cos(teta[2])*sin(teta[0])*sin(teta[3])) - 40*sin(teta[3])*(cos(teta[0])*sin(teta[2]) + cos(teta[1])*cos(teta[2])*sin(teta[0])) + 86*sin(teta[5])*(cos(teta[4])*(sin(teta[0])*sin(teta[1])*sin(teta[3]) + cos(teta[0])*cos(teta[3])*sin(teta[2]) + cos(teta[1])*cos(teta[2])*cos(teta[3])*sin(teta[0])) + sin(teta[4])*(cos(teta[0])*cos(teta[2]) - cos(teta[1])*sin(teta[0])*sin(teta[2]))) + 40*cos(teta[3])*sin(teta[0])*sin(teta[1]) - 75) - 8*(20*cos(teta[2])*cos(teta[3])*sin(teta[1]) - 20*cos(teta[1])*sin(teta[3]) - 43*cos(teta[1])*cos(teta[5])*sin(teta[3]) + 43*cos(teta[2])*cos(teta[3])*cos(teta[5])*sin(teta[1]) + 43*cos(teta[1])*cos(teta[3])*cos(teta[4])*sin(teta[5]) + 43*cos(teta[2])*cos(teta[4])*sin(teta[1])*sin(teta[3])*sin(teta[5]))*(20*cos(teta[1]) + 20*cos(teta[1])*cos(teta[3]) + 20*cos(teta[2])*sin(teta[1])*sin(teta[3]) + 43*cos(teta[1])*cos(teta[3])*cos(teta[5]) + 43*cos(teta[2])*cos(teta[5])*sin(teta[1])*sin(teta[3]) + 43*cos(teta[1])*cos(teta[4])*sin(teta[3])*sin(teta[5]) + 43*sin(teta[1])*sin(teta[2])*sin(teta[4])*sin(teta[5]) - 43*cos(teta[2])*cos(teta[3])*cos(teta[4])*sin(teta[1])*sin(teta[5]) + 7) + 2*(86*cos(teta[5])*(cos(teta[0])*sin(teta[1])*sin(teta[3]) - cos(teta[3])*sin(teta[0])*sin(teta[2]) + cos(teta[0])*cos(teta[1])*cos(teta[2])*cos(teta[3])) - 40*cos(teta[3])*(sin(teta[0])*sin(teta[2]) - cos(teta[0])*cos(teta[1])*cos(teta[2])) + 40*cos(teta[0])*sin(teta[1])*sin(teta[3]) - 86*cos(teta[4])*sin(teta[5])*(sin(teta[0])*sin(teta[2])*sin(teta[3]) + cos(teta[0])*cos(teta[3])*sin(teta[1]) - cos(teta[0])*cos(teta[1])*cos(teta[2])*sin(teta[3])))*(40*cos(teta[0])*sin(teta[1]) + 86*cos(teta[5])*(sin(teta[0])*sin(teta[2])*sin(teta[3]) + cos(teta[0])*cos(teta[3])*sin(teta[1]) - cos(teta[0])*cos(teta[1])*cos(teta[2])*sin(teta[3])) + 40*sin(teta[3])*(sin(teta[0])*sin(teta[2]) - cos(teta[0])*cos(teta[1])*cos(teta[2])) + 86*sin(teta[5])*(cos(teta[4])*(cos(teta[0])*sin(teta[1])*sin(teta[3]) - cos(teta[3])*sin(teta[0])*sin(teta[2]) + cos(teta[0])*cos(teta[1])*cos(teta[2])*cos(teta[3])) - sin(teta[4])*(cos(teta[2])*sin(teta[0]) + cos(teta[0])*cos(teta[1])*sin(teta[2]))) + 40*cos(teta[0])*cos(teta[3])*sin(teta[1]) - 75))/(200*sqrt( pow(40*sin(teta[0])*sin(teta[1]) - 86*cos(teta[5])*(cos(teta[0])*sin(teta[2])*sin(teta[3]) - cos(teta[3])*sin(teta[0])*sin(teta[1]) + cos(teta[1])*cos(teta[2])*sin(teta[0])*sin(teta[3])) - 40*sin(teta[3])*(cos(teta[0])*sin(teta[2]) + cos(teta[1])*cos(teta[2])*sin(teta[0])) + 86*sin(teta[5])*(cos(teta[4])*(sin(teta[0])*sin(teta[1])*sin(teta[3]) + cos(teta[0])*cos(teta[3])*sin(teta[2]) + cos(teta[1])*cos(teta[2])*cos(teta[3])*sin(teta[0])) + sin(teta[4])*(cos(teta[0])*cos(teta[2]) - cos(teta[1])*sin(teta[0])*sin(teta[2]))) + 40*cos(teta[3])*sin(teta[0])*sin(teta[1]) - 75,2) + 4* pow(20*cos(teta[1]) + 20*cos(teta[1])*cos(teta[3]) + 20*cos(teta[2])*sin(teta[1])*sin(teta[3]) + 43*cos(teta[1])*cos(teta[3])*cos(teta[5]) + 43*cos(teta[2])*cos(teta[5])*sin(teta[1])*sin(teta[3]) + 43*cos(teta[1])*cos(teta[4])*sin(teta[3])*sin(teta[5]) + 43*sin(teta[1])*sin(teta[2])*sin(teta[4])*sin(teta[5]) - 43*cos(teta[2])*cos(teta[3])*cos(teta[4])*sin(teta[1])*sin(teta[5]) + 7, 2) + pow(40*cos(teta[0])*sin(teta[1]) + 86*cos(teta[5])*(sin(teta[0])*sin(teta[2])*sin(teta[3]) + cos(teta[0])*cos(teta[3])*sin(teta[1]) - cos(teta[0])*cos(teta[1])*cos(teta[2])*sin(teta[3])) + 40*sin(teta[3])*(sin(teta[0])*sin(teta[2]) - cos(teta[0])*cos(teta[1])*cos(teta[2])) + 86*sin(teta[5])*(cos(teta[4])*(cos(teta[0])*sin(teta[1])*sin(teta[3]) - cos(teta[3])*sin(teta[0])*sin(teta[2]) + cos(teta[0])*cos(teta[1])*cos(teta[2])*cos(teta[3])) - sin(teta[4])*(cos(teta[2])*sin(teta[0]) + cos(teta[0])*cos(teta[1])*sin(teta[2]))) + 40*cos(teta[0])*cos(teta[3])*sin(teta[1]) - 75,2) ));


w_d[4]=-(172*sin(teta[5])*(sin(teta[4])*(cos(teta[0])*sin(teta[1])*sin(teta[3]) - cos(teta[3])*sin(teta[0])*sin(teta[2]) + cos(teta[0])*cos(teta[1])*cos(teta[2])*cos(teta[3])) + cos(teta[4])*(cos(teta[2])*sin(teta[0]) + cos(teta[0])*cos(teta[1])*sin(teta[2])))*(40*cos(teta[0])*sin(teta[1]) + 86*cos(teta[5])*(sin(teta[0])*sin(teta[2])*sin(teta[3]) + cos(teta[0])*cos(teta[3])*sin(teta[1]) - cos(teta[0])*cos(teta[1])*cos(teta[2])*sin(teta[3])) + 40*sin(teta[3])*(sin(teta[0])*sin(teta[2]) - cos(teta[0])*cos(teta[1])*cos(teta[2])) + 86*sin(teta[5])*(cos(teta[4])*(cos(teta[0])*sin(teta[1])*sin(teta[3]) - cos(teta[3])*sin(teta[0])*sin(teta[2]) + cos(teta[0])*cos(teta[1])*cos(teta[2])*cos(teta[3])) - sin(teta[4])*(cos(teta[2])*sin(teta[0]) + cos(teta[0])*cos(teta[1])*sin(teta[2]))) + 40*cos(teta[0])*cos(teta[3])*sin(teta[1]) - 75) + 172*sin(teta[5])*(sin(teta[4])*(sin(teta[0])*sin(teta[1])*sin(teta[3]) + cos(teta[0])*cos(teta[3])*sin(teta[2]) + cos(teta[1])*cos(teta[2])*cos(teta[3])*sin(teta[0])) - cos(teta[4])*(cos(teta[0])*cos(teta[2]) - cos(teta[1])*sin(teta[0])*sin(teta[2])))*(40*sin(teta[0])*sin(teta[1]) - 86*cos(teta[5])*(cos(teta[0])*sin(teta[2])*sin(teta[3]) - cos(teta[3])*sin(teta[0])*sin(teta[1]) + cos(teta[1])*cos(teta[2])*sin(teta[0])*sin(teta[3])) - 40*sin(teta[3])*(cos(teta[0])*sin(teta[2]) + cos(teta[1])*cos(teta[2])*sin(teta[0])) + 86*sin(teta[5])*(cos(teta[4])*(sin(teta[0])*sin(teta[1])*sin(teta[3]) + cos(teta[0])*cos(teta[3])*sin(teta[2]) + cos(teta[1])*cos(teta[2])*cos(teta[3])*sin(teta[0])) + sin(teta[4])*(cos(teta[0])*cos(teta[2]) - cos(teta[1])*sin(teta[0])*sin(teta[2]))) + 40*cos(teta[3])*sin(teta[0])*sin(teta[1]) - 75) - 344*sin(teta[5])*(cos(teta[4])*sin(teta[1])*sin(teta[2]) - cos(teta[1])*sin(teta[3])*sin(teta[4]) + cos(teta[2])*cos(teta[3])*sin(teta[1])*sin(teta[4]))*(20*cos(teta[1]) + 20*cos(teta[1])*cos(teta[3]) + 20*cos(teta[2])*sin(teta[1])*sin(teta[3]) + 43*cos(teta[1])*cos(teta[3])*cos(teta[5]) + 43*cos(teta[2])*cos(teta[5])*sin(teta[1])*sin(teta[3]) + 43*cos(teta[1])*cos(teta[4])*sin(teta[3])*sin(teta[5]) + 43*sin(teta[1])*sin(teta[2])*sin(teta[4])*sin(teta[5]) - 43*cos(teta[2])*cos(teta[3])*cos(teta[4])*sin(teta[1])*sin(teta[5]) + 7))/(200* sqrt(pow(40*sin(teta[0])*sin(teta[1]) - 86*cos(teta[5])*(cos(teta[0])*sin(teta[2])*sin(teta[3]) - cos(teta[3])*sin(teta[0])*sin(teta[1]) + cos(teta[1])*cos(teta[2])*sin(teta[0])*sin(teta[3])) - 40*sin(teta[3])*(cos(teta[0])*sin(teta[2]) + cos(teta[1])*cos(teta[2])*sin(teta[0])) + 86*sin(teta[5])*(cos(teta[4])*(sin(teta[0])*sin(teta[1])*sin(teta[3]) + cos(teta[0])*cos(teta[3])*sin(teta[2]) + cos(teta[1])*cos(teta[2])*cos(teta[3])*sin(teta[0])) + sin(teta[4])*(cos(teta[0])*cos(teta[2]) - cos(teta[1])*sin(teta[0])*sin(teta[2]))) + 40*cos(teta[3])*sin(teta[0])*sin(teta[1]) - 75,2) + 4* pow(20*cos(teta[1]) + 20*cos(teta[1])*cos(teta[3]) + 20*cos(teta[2])*sin(teta[1])*sin(teta[3]) + 43*cos(teta[1])*cos(teta[3])*cos(teta[5]) + 43*cos(teta[2])*cos(teta[5])*sin(teta[1])*sin(teta[3]) + 43*cos(teta[1])*cos(teta[4])*sin(teta[3])*sin(teta[5]) + 43*sin(teta[1])*sin(teta[2])*sin(teta[4])*sin(teta[5]) - 43*cos(teta[2])*cos(teta[3])*cos(teta[4])*sin(teta[1])*sin(teta[5]) + 7,2) + pow(40*cos(teta[0])*sin(teta[1]) + 86*cos(teta[5])*(sin(teta[0])*sin(teta[2])*sin(teta[3]) + cos(teta[0])*cos(teta[3])*sin(teta[1]) - cos(teta[0])*cos(teta[1])*cos(teta[2])*sin(teta[3])) + 40*sin(teta[3])*(sin(teta[0])*sin(teta[2]) - cos(teta[0])*cos(teta[1])*cos(teta[2])) + 86*sin(teta[5])*(cos(teta[4])*(cos(teta[0])*sin(teta[1])*sin(teta[3]) - cos(teta[3])*sin(teta[0])*sin(teta[2]) + cos(teta[0])*cos(teta[1])*cos(teta[2])*cos(teta[3])) - sin(teta[4])*(cos(teta[2])*sin(teta[0]) + cos(teta[0])*cos(teta[1])*sin(teta[2]))) + 40*cos(teta[0])*cos(teta[3])*sin(teta[1]) - 75,2)));



w_d[5] = -(8*(43*cos(teta[1])*cos(teta[3])*sin(teta[5]) - 43*cos(teta[1])*cos(teta[4])*cos(teta[5])*sin(teta[3]) + 43*cos(teta[2])*sin(teta[1])*sin(teta[3])*sin(teta[5]) - 43*cos(teta[5])*sin(teta[1])*sin(teta[2])*sin(teta[4]) + 43*cos(teta[2])*cos(teta[3])*cos(teta[4])*cos(teta[5])*sin(teta[1]))*(20*cos(teta[1]) + 20*cos(teta[1])*cos(teta[3]) + 20*cos(teta[2])*sin(teta[1])*sin(teta[3]) + 43*cos(teta[1])*cos(teta[3])*cos(teta[5]) + 43*cos(teta[2])*cos(teta[5])*sin(teta[1])*sin(teta[3]) + 43*cos(teta[1])*cos(teta[4])*sin(teta[3])*sin(teta[5]) + 43*sin(teta[1])*sin(teta[2])*sin(teta[4])*sin(teta[5]) - 43*cos(teta[2])*cos(teta[3])*cos(teta[4])*sin(teta[1])*sin(teta[5]) + 7) + 2*(86*sin(teta[5])*(sin(teta[0])*sin(teta[2])*sin(teta[3]) + cos(teta[0])*cos(teta[3])*sin(teta[1]) - cos(teta[0])*cos(teta[1])*cos(teta[2])*sin(teta[3])) - 86*cos(teta[5])*(cos(teta[4])*(cos(teta[0])*sin(teta[1])*sin(teta[3]) - cos(teta[3])*sin(teta[0])*sin(teta[2]) + cos(teta[0])*cos(teta[1])*cos(teta[2])*cos(teta[3])) - sin(teta[4])*(cos(teta[2])*sin(teta[0]) + cos(teta[0])*cos(teta[1])*sin(teta[2]))))*(40*cos(teta[0])*sin(teta[1]) + 86*cos(teta[5])*(sin(teta[0])*sin(teta[2])*sin(teta[3]) + cos(teta[0])*cos(teta[3])*sin(teta[1]) - cos(teta[0])*cos(teta[1])*cos(teta[2])*sin(teta[3])) + 40*sin(teta[3])*(sin(teta[0])*sin(teta[2]) - cos(teta[0])*cos(teta[1])*cos(teta[2])) + 86*sin(teta[5])*(cos(teta[4])*(cos(teta[0])*sin(teta[1])*sin(teta[3]) - cos(teta[3])*sin(teta[0])*sin(teta[2]) + cos(teta[0])*cos(teta[1])*cos(teta[2])*cos(teta[3])) - sin(teta[4])*(cos(teta[2])*sin(teta[0]) + cos(teta[0])*cos(teta[1])*sin(teta[2]))) + 40*cos(teta[0])*cos(teta[3])*sin(teta[1]) - 75) - 2*(86*sin(teta[5])*(cos(teta[0])*sin(teta[2])*sin(teta[3]) - cos(teta[3])*sin(teta[0])*sin(teta[1]) + cos(teta[1])*cos(teta[2])*sin(teta[0])*sin(teta[3])) + 86*cos(teta[5])*(cos(teta[4])*(sin(teta[0])*sin(teta[1])*sin(teta[3]) + cos(teta[0])*cos(teta[3])*sin(teta[2]) + cos(teta[1])*cos(teta[2])*cos(teta[3])*sin(teta[0])) + sin(teta[4])*(cos(teta[0])*cos(teta[2]) - cos(teta[1])*sin(teta[0])*sin(teta[2]))))*(40*sin(teta[0])*sin(teta[1]) - 86*cos(teta[5])*(cos(teta[0])*sin(teta[2])*sin(teta[3]) - cos(teta[3])*sin(teta[0])*sin(teta[1]) + cos(teta[1])*cos(teta[2])*sin(teta[0])*sin(teta[3])) - 40*sin(teta[3])*(cos(teta[0])*sin(teta[2]) + cos(teta[1])*cos(teta[2])*sin(teta[0])) + 86*sin(teta[5])*(cos(teta[4])*(sin(teta[0])*sin(teta[1])*sin(teta[3]) + cos(teta[0])*cos(teta[3])*sin(teta[2]) + cos(teta[1])*cos(teta[2])*cos(teta[3])*sin(teta[0])) + sin(teta[4])*(cos(teta[0])*cos(teta[2]) - cos(teta[1])*sin(teta[0])*sin(teta[2]))) + 40*cos(teta[3])*sin(teta[0])*sin(teta[1]) - 75))/(200*sqrt(pow(40*sin(teta[0])*sin(teta[1]) - 86*cos(teta[5])*(cos(teta[0])*sin(teta[2])*sin(teta[3]) - cos(teta[3])*sin(teta[0])*sin(teta[1]) + cos(teta[1])*cos(teta[2])*sin(teta[0])*sin(teta[3])) - 40*sin(teta[3])*(cos(teta[0])*sin(teta[2]) + cos(teta[1])*cos(teta[2])*sin(teta[0])) + 86*sin(teta[5])*(cos(teta[4])*(sin(teta[0])*sin(teta[1])*sin(teta[3]) + cos(teta[0])*cos(teta[3])*sin(teta[2]) + cos(teta[1])*cos(teta[2])*cos(teta[3])*sin(teta[0])) + sin(teta[4])*(cos(teta[0])*cos(teta[2]) - cos(teta[1])*sin(teta[0])*sin(teta[2]))) + 40*cos(teta[3])*sin(teta[0])*sin(teta[1]) - 75,2 )+ 4* pow(20*cos(teta[1]) + 20*cos(teta[1])*cos(teta[3]) + 20*cos(teta[2])*sin(teta[1])*sin(teta[3]) + 43*cos(teta[1])*cos(teta[3])*cos(teta[5]) + 43*cos(teta[2])*cos(teta[5])*sin(teta[1])*sin(teta[3]) + 43*cos(teta[1])*cos(teta[4])*sin(teta[3])*sin(teta[5]) + 43*sin(teta[1])*sin(teta[2])*sin(teta[4])*sin(teta[5]) - 43*cos(teta[2])*cos(teta[3])*cos(teta[4])*sin(teta[1])*sin(teta[5]) + 7,2 )+ pow(40*cos(teta[0])*sin(teta[1]) + 86*cos(teta[5])*(sin(teta[0])*sin(teta[2])*sin(teta[3]) + cos(teta[0])*cos(teta[3])*sin(teta[1]) - cos(teta[0])*cos(teta[1])*cos(teta[2])*sin(teta[3])) + 40*sin(teta[3])*(sin(teta[0])*sin(teta[2]) - cos(teta[0])*cos(teta[1])*cos(teta[2])) + 86*sin(teta[5])*(cos(teta[4])*(cos(teta[0])*sin(teta[1])*sin(teta[3]) - cos(teta[3])*sin(teta[0])*sin(teta[2]) + cos(teta[0])*cos(teta[1])*cos(teta[2])*cos(teta[3])) - sin(teta[4])*(cos(teta[2])*sin(teta[0]) + cos(teta[0])*cos(teta[1])*sin(teta[2]))) + 40*cos(teta[0])*cos(teta[3])*sin(teta[1]) - 75,2)));


w_d[6]=0;



}













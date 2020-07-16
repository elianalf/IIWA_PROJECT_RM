#include "ros/ros.h"
#include "boost/thread.hpp"
#include "sensor_msgs/JointState.h"
#include "geometry_msgs/Pose.h"
#include <std_msgs/Float64.h>
#include "nav_msgs/Path.h"
#include "KIN_LIB_IIWA.h"
#include "geometry_msgs/Twist.h"
#include <eigen3/Eigen/Geometry>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <iostream>
#include <fstream>
#include <tf2_geometry_msgs/tf2_geometry_msgs.h>
#include <stdio.h>
#include <string.h>



using namespace std;
   ofstream myfile;
   
class ADM_CONTROL {

   public:
      ADM_CONTROL();
      void run();
      void goto_work_pose();
      void joint_states_cb(sensor_msgs::JointState);
      void compute_dirkin();
      void ctrl_loop();
      
   private:
      ros::NodeHandle n;
     ros::Subscriber js_sub;
      ros::Publisher joint_pos_cmd[7];
      double q_cur[7];
      Eigen::Matrix<double,3,1> p_cur;
      Eigen::Matrix<double,4,1> o_cur;
      Eigen::Matrix4d T;
      std_msgs::Float64 q_cmd[7];
      int n_samples=2000; 
      Eigen::Matrix3d Rc;
		bool x_cur_available=false;
      bool new_torque = false;

};


//CONSTRUCTOR
ADM_CONTROL::ADM_CONTROL(){
   cout << "constructor"<<endl;
    joint_pos_cmd[0] = n.advertise<std_msgs::Float64>("/iiwa/joint1_position_controller/command", 0);
    joint_pos_cmd[1] = n.advertise<std_msgs::Float64>("/iiwa/joint2_position_controller/command", 0);
    joint_pos_cmd[2] = n.advertise<std_msgs::Float64>("/iiwa/joint3_position_controller/command", 0);
    joint_pos_cmd[3] = n.advertise<std_msgs::Float64>("/iiwa/joint4_position_controller/command", 0);
    joint_pos_cmd[4] = n.advertise<std_msgs::Float64>("/iiwa/joint5_position_controller/command", 0);
    joint_pos_cmd[5] = n.advertise<std_msgs::Float64>("/iiwa/joint6_position_controller/command", 0);
    joint_pos_cmd[6] = n.advertise<std_msgs::Float64>("/iiwa/joint7_position_controller/command", 0); 
    js_sub = n.subscribe("/iiwa/joint_states", 0, &ADM_CONTROL::joint_states_cb, this);
    p_cur.setZero();
	 goto_work_pose();
	cout<<"Robot ready to work"<<endl;
	new_torque=true;
	
}

//GO TO START POSE
void ADM_CONTROL::goto_work_pose() {
	double q_data[7];
	ros::Rate r(50);
	
   q_data[0] = 0;
	q_data[1] = 0.92; //x2
	q_data[2] = 0.0; //z3
	q_data[3] = -1.42; //x4
	q_data[4] = 0.0;//z5
	q_data[5] = 0.81;//x6
	q_data[6] = 0.0;//z7
	 
	for (int j=0; j<10;j++){
      for(int i=0; i<7; i++){
          q_cmd[i].data = q_data[i];
          joint_pos_cmd[i].publish(q_cmd[i]);
      }
      r.sleep();
   }
   sleep(1);
}


//GET CURRENT JOINT STATES 
void ADM_CONTROL::joint_states_cb(sensor_msgs::JointState js){
   for(int i=0; i<7; i++ ){
		q_cur[i] = js.position[i];
      //cout<<"Current joint states: "<< q_cur[i];
   }
}


Eigen::Matrix<double,3,3> skew(Eigen::Matrix<double,3,1> v)
{
	Eigen::Matrix<double,3,3> S;
	S << 0.0,		-v[2],	 v[1],		//Skew-symmetric matrix
						 v[2],	 0.0,	-v[0],
						-v[1],	 v[0],  0.0;
	return S;
}


Eigen::Matrix<double,3,1> QuatErr(Eigen::Matrix<double,4,1> Qd, Eigen::Matrix<double,4,1> Qe)
{
	Eigen::Matrix<double,4,4> Rd;
	Eigen::Matrix<double,4,4> Re;
	Eigen::Matrix<double,3,1> eps_d; 
	Eigen::Matrix<double,3,1> eps_e;
	Eigen::Matrix<double,3,1> eo;
	Eigen::Matrix<double,3,3> S;
	eps_d <<Qd[1], Qd[2], Qd[3];
	eps_e <<Qe[1], Qe[2], Qe[3];
	S = skew(eps_d);
	eo = Qe[0] * eps_d - Qd[0] * eps_e - S * eps_e;
	return eo;
	
}

Eigen::Matrix<double,4,1> quat_prop_function( Eigen::Matrix<double,3,1> w, Eigen::Matrix<double,4,1> quat_p ) {
  	double eta = quat_p[0];
  	Eigen::Matrix<double,1,3> eps;
  	eps << quat_p[1],  quat_p[2],  quat_p[3] ;
  	Eigen::Matrix<double,3,3> I;
  	I.setIdentity();
  	double deta = -0.5*eps*w;
  	Eigen::Matrix<double,3,3> S = skew( eps );
  	Eigen::Matrix<double,3,1> deps; 
  	deps = 0.5*( eta*I-S)*w;
   Eigen::Matrix<double,4,1> Q_P;
   Q_P << deta, deps[0], deps[1], deps[2] ;
  	return Q_P;
} 



//DIRECT KINEMATICS
void ADM_CONTROL::compute_dirkin() {
   KIN_LIB_IIWA cl;
	Eigen::Matrix4d T;
	
	while( !new_torque ) usleep(100);

	ros::Rate r(801);
	while(ros::ok()) {
	   
	   cl.direct_kin(q_cur, T);
	   
	   Rc = T.block<3,3>(0,0);
	   Eigen::Quaterniond quat(Rc);
	   quat.normalize();
	   
	   o_cur[0] = quat.w();
	   o_cur[1] = quat.x();
	   o_cur[2] = quat.y();
	   o_cur[3] = quat.z();
	   
	   p_cur = T.block<3,1>(0,3);
	   
	   //cout<<"o_cur "<<o_cur<<endl;
	 
	   r.sleep();
	   
	  x_cur_available=true;
	}
}


//ADMITTANCE CONTROL AND CLIK 
void ADM_CONTROL::ctrl_loop(){
  // cout<<"wait"<<endl;
   double q_data[7];
   
   double dt=0.00125;
   while( !x_cur_available) {usleep(100);}
  // cout<<"c_loop"<<endl;
  Eigen::Matrix< double, 6, 1> m ;
   Eigen::Matrix< double, 6, 1> k ;
   Eigen::Matrix< double, 6, 1> kp ;
   Eigen::Matrix< double, 6, 1> kd ;
    k << 100, 100, 100, 1, 1, 1;
    m << 10, 10, 10, 1, 1, 1;
    kp <<0, 0, 0, 0, 0, 0; 
    //kp <<700, 700, 700, 500, 500, 500;
    kd << 500,500, 500, 500, 500, 500;
   Eigen::Matrix< double, 6, 6> Kp = kp.asDiagonal();
   Eigen::Matrix< double, 6, 6> K = k.asDiagonal();
   Eigen::Matrix< double, 6, 6> Kd = kd.asDiagonal();
   Eigen::Matrix< double, 6, 6> M = m.asDiagonal();
   Eigen::Matrix<double,6,7> J;
   KIN_LIB_IIWA cl;
   Eigen::Matrix<double,6,1> he;
   Eigen::Matrix<double,7,1> dq;
   Eigen::Matrix<double,7,6> Jpinv;
   Eigen::Matrix<double,6,1> e;
   Eigen::Matrix<double,3,1> eo;
   Eigen::Matrix<double,3,1> x_i;
   Eigen::Matrix<double,3,1> p_z;
   Eigen::Matrix<double,6,1> z_tilde;
   Eigen::Matrix<double,6,1> z_d;
   Eigen::Matrix<double,6,1> z_dd;
   Eigen::Matrix<double,2000,1> sx;
   Eigen::Matrix<double,2000,1> sy;
   Eigen::Matrix<double,2000,1> sz;
    double roll[2000];
    double pitch[2000];
    double yaw[2000];
    Eigen::Vector3d euler_angles;  
    Eigen::Matrix<double,4,1> o_i;
    Eigen::Matrix<double,4,1> o_z;
    Eigen::Matrix<double,4,1> quat_prop;
    Eigen::Matrix<double,3,1> z_tilde_o;
    z_tilde_o.setZero();
   o_i.setZero();
   p_z.setZero();
   z_d.setZero();
   z_dd.setZero();
   x_i.setZero();
   dq.setZero();
   e.setZero();
   he.setZero();
   sx.setZero();
   sy.setZero();
   sz.setZero();
   int index=0;
   
   ros::Rate c_rate(800);

   o_i = o_cur;
   x_i = p_cur;
   p_z = p_cur;
   o_z=o_i;
  while(index<n_samples){
      he[0] = 0;
      he[1]= 0;
      he[2]= 0;
     if((index>50)&&(index<850))
         { 
            he[2]= 15;// 2 = z; 4 rot intorno y
         }
    /* if((index>50)&&(index<850))
         { double _arg ;
            _arg=2*3.14*(index%400)/400;
            he[2] = 15*sin(_arg);
         }*/
    cout<<"HE: "<<he<<endl;
    
     z_tilde_o=QuatErr(o_i,o_z);    
     z_tilde.topRows<3>()=(x_i-p_z);
     z_tilde.bottomRows<3>()=z_tilde_o;
     
     z_dd = M.inverse()*(he - Kd*z_d + Kp*z_tilde);
      //cout<<"z_dd "<<z_dd<<endl;
      
     z_d = z_dd*dt + z_d;
     //cout<<"z_d "<<z_d<<endl;
     
     quat_prop = quat_prop_function( z_d.bottomRows<3>(), o_z); 
     o_z = quat_prop*dt + o_z;
     o_z.normalize();

     p_z = z_d.topRows<3>() *dt + p_z;
     
      cl.compute_jacobian(q_cur,J);
      
	   e.topRows<3>() = p_z - p_cur;  
	   eo=QuatErr(o_z,o_cur);
      e.bottomRows<3>()=eo;
	   //cout<<"e:"<<endl;
      //cout<<e<<endl;


     Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cqr(J);
      Jpinv= cqr.pseudoInverse();
      dq= Jpinv*(z_d + K * e);
     // cout<<"q_dot:"<<dq<<endl;
      
      for(int j=0;j<7;j++){
         q_data[j] = dq[j]*dt + q_cur[j];
	   }
	   cout<<"PUBLISH COMMAND "<<index<<endl;
      for(int i=0; i<7; i++){
          q_cmd[i].data = q_data[i];
          joint_pos_cmd[i].publish(q_cmd[i]);
          //cout<<q_data[i] <<" ";
      }
     //cout<<endl;
     
     // **Vectors to save the results: **
     sx[index]=p_cur[0]-x_i[0];
     sy[index]=p_cur[1]-x_i[1];
     sz[index]=p_cur[2]-x_i[2];
      //cout<<"s: "<<sx[index]<<" "<<sy[index]<<" "<<sz[index]<<endl;
      euler_angles = Rc.eulerAngles(0, 1, 2);
      roll[index]=euler_angles[0];
      pitch[index]=euler_angles[1];
      yaw[index]=euler_angles[2];
      c_rate.sleep();
       
      index++;
   
   }

   
 myfile.open ("adm_spost5.txt");
if (!myfile.is_open()) {cout<<"*************ERROR****************";}
  myfile << "[ "; 
     for(int i=0;i < n_samples;i++){
        myfile << sx[i];
        myfile << "; ";
   }
   myfile << " ] \n"; 
   myfile << "[ "; 
     for(int i=0;i <  n_samples;i++){
        myfile <<sy[i];
        myfile << "; ";
   }
   myfile << " ] \n"; 
   myfile << "[ "; 
     for(int i=0;i < n_samples;i++){
        myfile <<sz[i];
        myfile << "; ";
   }
   myfile << " ] \n"; 
   myfile.close();
    myfile.open ("adm_orient5.txt");
   if (!myfile.is_open()) {cout<<"*************ERROR****************";}
  myfile << "[ "; 
     for(int i=0;i < n_samples;i++){
        myfile << roll[i];
        myfile << "; ";
   }
   myfile << " ] \n"; 
   myfile << "[ "; 
     for(int i=0;i <  n_samples;i++){
        myfile <<pitch[i];
        myfile << "; ";
   }
   myfile << " ] \n"; 
   myfile << "[ "; 
     for(int i=0;i < n_samples;i++){
        myfile <<yaw[i];
        myfile << "; ";
   }
   myfile << " ] \n";
    myfile.close();
   cout<<"END"<<endl;
   
}

void ADM_CONTROL::run(){
   cout<< "run"<<endl;
   boost::thread get_dirkin_t( &ADM_CONTROL::compute_dirkin, this);
   boost::thread ctrl_loop_t ( &ADM_CONTROL::ctrl_loop, this);
   ros::spin();

}


int main(int argc, char** argv) {

	ros::init(argc, argv, "iiwa_admittance_control");
	ADM_CONTROL adc;
	adc.run();

	return 0;
}

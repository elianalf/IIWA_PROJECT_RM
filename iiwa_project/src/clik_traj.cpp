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

using namespace std;
   ofstream myfile;
 
class CLIK_ {

   public:
      CLIK_();
      void run();
      void goto_work_pose();
      void get_path();
      void joint_states_cb(sensor_msgs::JointState);
      void compute_dirkin();
      void ctrl_loop();
     void  get_des_velocity();
   private:
      ros::NodeHandle n;
      ros::Subscriber js_sub;
      ros::Publisher joint_pos_cmd[7];
      double q_cur[7];
      Eigen::Matrix3d Rc;
      Eigen::Matrix<double,3,1> p_cur;
      Eigen::Matrix<double,4,1> o_cur;
      Eigen::Matrix4d T;
      std_msgs::Float64 q_cmd[7];
      nav_msgs::Path ref_path;
      int _path_size=2100; 
      geometry_msgs::Twist des_velocity[2100];
	bool x_cur_available=false;
      bool new_path = false;
      
};


//CLASS CONSTRUCTOR
CLIK_::CLIK_(){ 
	cout<<"Begin clik"<<endl;
   joint_pos_cmd[0] = n.advertise<std_msgs::Float64>("/iiwa/joint1_position_controller/command", 0);
    joint_pos_cmd[1] = n.advertise<std_msgs::Float64>("/iiwa/joint2_position_controller/command", 0);
    joint_pos_cmd[2] = n.advertise<std_msgs::Float64>("/iiwa/joint3_position_controller/command", 0);
    joint_pos_cmd[3] = n.advertise<std_msgs::Float64>("/iiwa/joint4_position_controller/command", 0);
    joint_pos_cmd[4] = n.advertise<std_msgs::Float64>("/iiwa/joint5_position_controller/command", 0);
    joint_pos_cmd[5] = n.advertise<std_msgs::Float64>("/iiwa/joint6_position_controller/command", 0);
    joint_pos_cmd[6] = n.advertise<std_msgs::Float64>("/iiwa/joint7_position_controller/command", 0); 
    js_sub = n.subscribe("/iiwa/joint_states", 0, &CLIK_::joint_states_cb, this);
    
	goto_work_pose();
	get_des_velocity();
   get_path();
	cout<<"Robot ready to work"<<endl;
	usleep(100000);
   new_path = true;
   
}


Eigen::Matrix<double,3,3> skew(Eigen::Matrix<double,3,1> v)
{
	Eigen::Matrix<double,3,3> S;
	S << 0.0,		-v[2],	 v[1],		//Skew-symmetric matrix
						 v[2],		 0.0,	-v[0],
						-v[1],		 v[0], 	 0.0;
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


//GO TO START POSE 
void CLIK_::goto_work_pose() {
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
   usleep(1000000);
}


//GET CURRENT JOINT STATES
void CLIK_::joint_states_cb(sensor_msgs::JointState js){
   for(int i=0; i<7; i++ ){
		q_cur[i] = js.position[i];
      //cout<<"Current joint states: "<< q_cur[i];
   }
}


//GET VELOCITY FROM FILE
void CLIK_::get_des_velocity(){
     string line2;
  ifstream myfile2 ("ped_800.txt");
  int vindex=0;
  int indexv=0;

  if (myfile2.is_open())
  {
    while ( getline (myfile2,line2) )
    { 
      if(line2.at(0) == 'y'){
        vindex=1;
         indexv=0;}
      else if(line2.at(0) == 'z'){
         vindex=2;
         indexv=0;}
      else if(vindex==0){
         std::string::size_type sz;
        // cout << line2 << '\n';
         des_velocity[indexv].linear.x = std::stod(line2, &sz);
         // cout<<"x:"<<des_velocity[indexv].linear.x<<" i: "<<indexv<<endl;
          usleep(10);
          des_velocity[indexv].angular.x =0;
        des_velocity[indexv].angular.y =0;
       des_velocity[indexv].angular.z =0;
          indexv++;
      }
      else if(vindex==1){
         std::string::size_type sz;
       //  cout << line2 << '\n';
         des_velocity[indexv].linear.y = std::stod(line2, &sz);
         //cout<<"y:"<<des_velocity[indexv].linear.y<<" i: "<<indexv<<endl;
          usleep(10);
          indexv++;
      }
      else if(vindex==2){
      std::string::size_type sz;
         //cout << line2 << '\n';
         des_velocity[indexv].linear.z = std::stod(line2, &sz);
       //  cout<<"z:"<<des_velocity[indexv].linear.z<<" i: "<<indexv<<endl;
          usleep(10);
          indexv++;
      }
    }
    myfile2.close();
  }
    else cout << "Unable to open file"; 
     
     // cout<< "vel: "<< ve.linear.x<<" "<<ve.linear.y<< " "<<ve.linear.z<<" "<<index_vel<<endl; 
     //cout<< "vel: "<< des_velocity[index_vel].linear.x<<" "<<des_velocity[index_vel].linear.y<< " "<<des_velocity[index_vel].linear.z<<" "<<index_vel<<endl; 
  
}


//GET PATH FROM FILE 
void CLIK_::get_path(){
  string line;
  ifstream myfile ("pe_800.txt");
  double pro;
  int x_y_z=0;
  int i_myfile=0;
  
  geometry_msgs::PoseStamped my_path[2100];
  if (myfile.is_open())
  { cout<<"file opened"<<endl;
    while ( getline (myfile,line) )
    { 
      if(line.at(0) == 'y'){
         x_y_z=1;
         i_myfile=0;}
      else if(line.at(0) == 'z'){
         x_y_z=2;
         i_myfile=0;}
      else if(x_y_z==0){
         std::string::size_type sz;
        // cout << line << '\n';
         my_path[i_myfile].pose.position.x = std::stod(line, &sz);
          //cout<<"x:"<<my_path[i_myfile].pose.position.x<<" i: "<<i_myfile<<endl;
          usleep(10);
          my_path[i_myfile].pose.orientation.x =-0.219584460492;
         my_path[i_myfile].pose.orientation.y=0.97124257295;
         my_path[i_myfile].pose.orientation.z=-0.090106047927;
        my_path[i_myfile].pose.orientation.w=0.018744756398;
          i_myfile++;
      }
      else if(x_y_z==1){
         std::string::size_type sz;
        // cout << line << '\n';
         my_path[i_myfile].pose.position.y = std::stod(line, &sz);
         //cout<<"y:"<<my_path[i_myfile].pose.position.y<<" i: "<<i_myfile<<endl;
          usleep(10);
          i_myfile++;
      }
      else if(x_y_z==2){
      std::string::size_type sz;
         //cout << line << '\n';
         my_path[i_myfile].pose.position.z = std::stod(line, &sz);
         //cout<<"z:"<<my_path[i_myfile].pose.position.z<<" i: "<<i_myfile<<endl;
          usleep(10);
          i_myfile++;
      }
     
    }
    myfile.close();
  }

  else cout << "Unable to open file"; 
  
   for(int i=0;i < _path_size;i++){
     ref_path.poses.push_back(my_path[i]);
     //  cout<< "Point: "<< ref_path.poses[i].pose.position.x<<" "<<ref_path.poses[i].pose.position.y<< " "<<ref_path.poses[i].pose.position.z<<endl;
  }
  
      cout<<"path saved"<<endl;
    
}


//DIRECT KINEMATICS
void CLIK_::compute_dirkin() {
   KIN_LIB_IIWA cl;
	Eigen::Matrix4d T;
	
  
	while( !new_path ) usleep(100);
	cout<<"DIR KIN"<<endl;
	
	ros::Rate r(801);
	while(ros::ok()) {
	   cl.direct_kin(q_cur, T);
	 
	   Rc = T.block<3,3>(0,0);
	   Eigen::Quaterniond quat(Rc);
	   o_cur[0] = quat.w();
	    o_cur[1]=quat.x();
	    o_cur[2] =quat.y() ;
	    o_cur[3]= quat.z();
	   p_cur = T.block<3,1>(0,3);
	   //cout<<"p DIR: "<<p_cur[0]<<" "<<p_cur[1]<<" "<<p_cur[2]<<endl;
	   r.sleep();
	   x_cur_available=true;
	}
}


//CONTROL LOOP: CLIK
void CLIK_::ctrl_loop(){
   //cout<<"wait"<<endl;
   double q_data[7];
   double dt=0.00125;
   while( !x_cur_available) {usleep(100);}
   cout<<"c_loop"<<endl;
   usleep(100000);
   Eigen::Matrix< double, 6, 1> k ;
    k << 1450, 1450, 1450,1, 1,1;
    Eigen::Matrix< double, 6, 6> Kp = k.asDiagonal();
   double q_current[7];
   ros::Rate c_rate(800);
   Eigen::Matrix<double,6,7> J;
   KIN_LIB_IIWA cl;
   Eigen::Matrix<double,7,1> q0;
   Eigen::Matrix<double,6,1> x_dot;
   Eigen::Matrix<double,7,1> dq;
   Eigen::Matrix<double,7,6> Jpinv;
   Eigen::Matrix<double,6,1> e;
   Eigen::Matrix<double,3,1> eo;
   Eigen::Matrix3d Rd;
   double w_dq[7];
   Eigen::Matrix<double,2100,1> ex;
   Eigen::Matrix<double,2100,1> ey;
   Eigen::Matrix<double,2100,1> ez;
   Eigen::Matrix<double,4,1> o;
   Eigen::Matrix<double,7,7> I;
   I.setIdentity();
   o.setZero();
   e.setZero();
   eo.setZero();
   x_dot.setZero();
        
   for(int index=0; index< _path_size; index++){
      //cout<<"index: "<<index<<endl;
      x_dot[0]=des_velocity[index].linear.x ;
      x_dot[1]=des_velocity[index].linear.y;
      x_dot[2]=des_velocity[index].linear.z ;
      x_dot[3]=des_velocity[index].angular.x ;
      x_dot[4]=des_velocity[index].angular.y ;
      x_dot[5]=des_velocity[index].angular.z;
     for (int j=0; j<7;j++)
         q_current[j]=q_cur[j];
         
       cl.compute_jacobian(q_current,J);
	   //cout<<"Jacobian:"<<endl;
      //cout<<J<<endl;
      //cout<<"p CONTR: "<<p_cur[0]<<" "<<p_cur[1]<<" "<<p_cur[2]<<endl;
      //cout<< "vel: "<< des_velocity[index].linear.x<<" "<<des_velocity[index].linear.y<< " "<<des_velocity[index].linear.z<<endl;
      //cout<<"v angl: "<< des_velocity[index].angular.x <<des_velocity[index].angular.y<<des_velocity[index].angular.z<<endl;
      //cout<< "POINT: "<< ref_path.poses[index].pose.position.x<<" "<<ref_path.poses[index].pose.position.y<< " "<<ref_path.poses[index].pose.position.z<<endl;
      
      //**SECOND TASK** 
      //cl.compute_w_dq(q_current, w_dq);
     // for (int j=0; j<7;j++)
      //   q0[j]=1*w_dq[j];
      
      e[0]=ref_path.poses[index].pose.position.x - p_cur[0];
      e[1]=ref_path.poses[index].pose.position.y - p_cur[1];
      e[2]=ref_path.poses[index].pose.position.z - p_cur[2];
      ex[index]=e[0];
      ey[index]=e[1];
      ez[index]=e[2];
      Rd=Rc;
     Eigen::Quaterniond quat(Rd);
     quat.normalize();
     o[0] = quat.w();
     o[1] = quat.x();
     o[2] = quat.y();
     o[3] = quat.z();
     eo=QuatErr(o,o_cur);
     e.bottomRows<3>()=eo;
     //cout<<"e o "<<eo<<endl;
    //cout<<"e: "<<e[0]<<" "<<e[1]<<" "<<e[2]<<endl;
    
      Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cqr(J);
      Jpinv= cqr.pseudoInverse();
      //dq= Jpinv*(x_dot + Kp * e) + ( I - Jpinv*J)*q0;
      dq= Jpinv*(x_dot + Kp * e);
      //cout<<"q_dot:"<<dq<<endl;
      
      for(int j=0;j<7;j++){
         q_data[j] = dq[j]*dt + q_current[j];
	   }
	   cout<<"PUBLISH COMMAND"<<endl;
      for(int i=0; i<7; i++){
          q_cmd[i].data = q_data[i];
          joint_pos_cmd[i].publish(q_cmd[i]);
         // cout<<q_data[i] <<" ";
      }
      //cout<<endl;
      c_rate.sleep();
 
   }
   myfile.open ("900hz_1500.txt");
   if (!myfile.is_open()) {cout<<"*************ERROR****************";}
    myfile << "[ "; 
     for(int i=0;i < _path_size;i++){
        myfile << ex[i];
        myfile << "; ";
   }
   myfile << " ] \n"; 
   myfile << "[ "; 
     for(int i=0;i <  _path_size;i++){
        myfile <<ey[i];
        myfile << "; ";
   }
   myfile << " ] \n"; 
   myfile << "[ "; 
     for(int i=0;i < _path_size;i++){
        myfile <<ez[i];
        myfile << "; ";
   }
   myfile << " ] \n"; 
    myfile.close();
    cout<<"END..."<<endl;
}


void CLIK_::run(){
   
   
   boost::thread ctrl_loop_t ( &CLIK_::ctrl_loop, this);
   boost::thread get_dirkin_t( &CLIK_::compute_dirkin, this);
   ros::spin();
}


int main(int argc, char** argv){
	ros::init(argc, argv, "iiwa_CLIK_");
	CLIK_ cl;
	cl.run();
	return 0;
}

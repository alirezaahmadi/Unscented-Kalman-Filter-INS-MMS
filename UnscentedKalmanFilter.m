%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UnscentedKalman Fliter  on IMU+GNSS Recorded Data .    %
% by: Alireza Ahmadi                                     %
% University of Bonn- MSc Robotics & Geodetic Engineering%
% Alireza.Ahmadi@uni-bonn.de                             %
% AhmadiAlireza.webs.com                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;
%%***************definition of parameters************************
 a = 6378137.0;        %m Earth's ellipsoid radius at wquator
 b = 6356752.3142 ;    %m Earth's ellipsoid radius at poles
 ecc = 0.0818191908426;  %- Earth's ellipsoid eccentricity
 w_i_e = 7.2921150*10^-5;   %rad/s Earth's rotational speed
 mu = 3.986004418*10^14;  %m^3/s^2 Geocentric gravitional constant
 f = 1/298.257223563;  %- Earth's ellipsoid flattening
 omega_ib_b = zeros(3,3);
 g0 = 0;
 R2D = 180/pi;
 D2R = pi/180;
%%***************************************************************
%filename0  = 'IMAR0000.mat';
%filename1  = 'IMAR0001.mat';
filename2  = 'IMAR0002.mat';
%filename3  = 'IMAR0003.mat';

load(filename2,'-regexp','imu');
load(filename2,'-regexp','gps');

%filename00  = 'UTM_IMAR0000.mat';
%filename01  = 'UTM_IMAR0001.mat';
filename02  = 'UTM_IMAR0002.mat';
%filename03  = 'UTM_IMAR0003.mat';

load(filename02,'-regexp','UTM');

% Initialization of State vector - first position is gotten from GPS
% output and firs heading come from compass data - others accelerations and
% velocities as the vehicle soppused to be stationary are zero ! (t = 0)
x_0 = UTM(1,1);
y_0 = UTM(1,2);
phi_0 = imu.rpy_ned(1,3);
phi_D_0 = 0;
v_0 = 0;
a_0 = 0;

% Perior state vector
x_k = [x_0,y_0,phi_0,phi_D_0,v_0,a_0]';

% State vecotr parameters in time (t+1)
x_i = 0;
y_i = 0;
phi_i = 0;
phi_D_i = 0;
v_i = 0;
a_i = 0;
% Posterior state vector
X_K = [x_i,y_i,phi_i,phi_D_i,v_0,a_i]';

% Time Step of proceeding the loop
Delta_t=0.001;

% Uncertainity Variance of sesnsors based on Datasheet of sensors !!
Variance_phi_D = 10^-12;
Varince_a = 10^-7;

% Estimated Accuracy of GPS data for recorded positions and heading
Variance_x = 10^-4;
Variance_y = 10^-4;
Variance_phi = 10^-5;

% Initial state covariance Matix

P_k = [Variance_x,0,0,0,0,0;
       0,Variance_y,0,0,0,0;
       0,0,Variance_phi,0,0,0;
       0,0,0,1,0,0;
       0,0,0,0,1,0;
       0,0,0,0,0,1];
%P_k = eye(6);
       
% some hardcoded counters to instantiate the loop
A = size(UTM);
Res = A(1,1);
B = size(imu);
cnt = 1;
cnt_b = 1;
Interval =1;
P_Interval =10;

% Observation Variance factors
Varince_x_gps = 0.01;
Varince_y_gps = 0.01;
Variance_Acc  = 0.001;
Variance_gyro = 5*10^-8;

% Initial Plot of GNSS recorded Positions (mm-accuracy)
% we print one position per each 15 recorded data to speedup the process as
% reference of comparion to sensor fusion and filters performance
for cnt=1:15:length(UTM)
P1 =plot(UTM(cnt,1),UTM(cnt,2),'r.');
hold on
end


cnt=Interval+1;
while  cnt<=Res
cnt;

%*****************Prediction*****************************
% State Vector
% f_x = [x_0 + v_0*Delta_t*cos(phi_0);
%       y_0 + v_0*Delta_t*sin(phi_0);
%       phi_0 + phi_D_0*Delta_t;
%       phi_D_0;
%       v_0+a_0*Delta_t;
%       a_0];
f_x=@(x)[x(1) + x(5) * cos(x(3) * Delta_t);
         x(2) + x(5) * sin(x(3) * Delta_t);
         x(3) + x(4) * Delta_t;
         x(4);
         x(5) + x(6) * Delta_t;
         x(6)];
% Command noise matrix
% Q = [Variance_phi_D,0,0,0,0,0;
%      0,     Varince_a,0,0,0,0;
%      0,0,0,0,0,0;
%      0,0,0,0,0,0;
%      0,0,0,0,0,0;
%      0,0,0,0,0,0;];
Q = [0,0,0,0,0,0;
     0,0,0,0,0,0;
     0,0,0,0,0,0;
     0,0,0,Variance_phi_D,0,0;
     0,0,0,0,0,0;
     0,0,0,0,0,Varince_a;];
   
% Observation Vector
z_k = [UTM(cnt,1);
       UTM(cnt,2);
       -imu.acc_ib_b(cnt_b,1);
       imu.omg_ib_b(cnt_b,3)];
   
% State to measurement mapping matrix
h_z = @(x)[x(1);
           x(2);
          -x(3);
           x(4)];

% Observation Noise matrix - all sensors assumed to be uncorrelated 
R = [Varince_x_gps,0,0,0;
     0,Varince_y_gps,0,0;
     0,0,Variance_Acc, 0;
     0,0,0,Variance_gyro];
       
% Unscented Kalman filter function
[X_K,P_k] = unscentedKF(f_x,x_k,P_k,h_z,z_k,Q,R);

%C_k = P_k * f_x * P_k^-1;
%x_kn = X_K + C_k(x_K1n - x_K1k);
%p_kn = P_K + C_k(p_K1n - p_K1k);

hold on
P2=plot(X_K(1,1),X_K(2,1),'b.');
  
cnt = cnt+Interval;
cnt_b = cnt_b + P_Interval;
x_k = X_K;
%Time fraction 
Delta_t = (imu.imu_time(cnt_b,1) - imu.imu_time(cnt_b-P_Interval,1));

end
xlabel('UTM-East');
ylabel('UTM-North');
legend([P1,P2],'UTM','UKF');
grid on;












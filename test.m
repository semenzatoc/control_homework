%% Reset
clearvars; clear memory;  clc; warning off
close all;
addpath(genpath('src'));
load('continuousModel.mat')
load('operatingPoint.mat')

%% Load operating point
x_eq=operatingPoint.getstatestruct.signals.values;
x_eq(1:4)=operatingPoint.States(1,1).x;
x_eq(5)=operatingPoint.States(2,1).x;
x_eq(6)=operatingPoint.States(3,1).x;
x_eq(7:10)=operatingPoint.States(4,1).x;
bis_eq=50;
rass_eq=-4;
propofol_eq=operatingPoint.getinputstruct.signals(1).values;
remifentanil_eq=operatingPoint.getinputstruct.signals(2).values;
op.propofol_eq=propofol_eq;
op.remifentanil_eq=remifentanil_eq;
op.bis_eq=bis_eq;
op.rass_eq=rass_eq;

%%    PARAMETERS FOR THE SPECIFIC SIMULATION
info.Tsim = 6000; % [s] simulation time
info.Ts=20;% [s] assumed to be equal for all signals (BIS, RASS and NMB)

wantNoisyBIS = 1;               % 0= I do NOT want noise on BIS, 1= I want noise on BIS
% wantNoisyNBM = 1;               % 0= I do NOT want noise on NBM, 1= I want noise on NBM

info.T_induction=0;
%BIS
info.initial_BIS=100;
info.desired_BIS = 50;
info.upperLimit_BIS= 60;
info.lowerLimit_BIS= 40;

%RASS
info.initial_RASS=0;
info.desired_RASS = -4;
info.upperLimit_RASS= -4;
info.lowerLimit_RASS= -5;

%% Patient database
subj = 0; % patient of interest MAX 24, 0 for the average subject
[patient] = patient_parameters(subj);

%% Drug parameters and hemodynamics information
[propofol, remifentanil, RASS] = drugs_parameters(patient);

%% Simulink simulation
% open loop analysis
% info.PropofolDose = 0.25;    % [0 5] 
% info.RemifentanilDose = 1.5; % [0 25]
info.PropofolDose = propofol_eq;
info.RemifentanilDose = remifentanil_eq;

%% Discretization 
discreteModel = c2d(continuousModel, info.Ts);

%% TASK 2 - MPC CONTROLLER
u_eq = [propofol_eq; remifentanil_eq];

discreteMat.A = discreteModel.A;
discreteMat.B = discreteModel.B;
discreteMat.M = discreteModel.B*0;
discreteMat.C = discreteModel.C;
discreteMat.D = discreteModel.D;

q=10*eye(2); 
r= 0.01*eye(2);
PH=200;
F= -[1 0; 0 1];
f= u_eq;

condDiscMat=...
BuildCondendsedMPCmatrices(discreteMat,r,q,F,f,PH);

%% KALMAN FILTER
Q_KF = diag((0.1*x_eq + min(x_eq)).^2);
R_KF = 0.001*eye(2);

% %Setting the new parameter for the MPC
q_Kalman=1*eye(2); 
r_Kalman= 0.01*eye(2);

kalmDiscMat=...
BuildCondendsedMPCmatrices(discreteMat,r_Kalman,q_Kalman,F,f,PH);

sim('PatientSIMULATOR_MIMOnoNMB')

% FIGURE: MPC Controller - Accessible state
ax=[];
figure('color', 'w'); 
sgtitle('MPC Control with accessible state')
ax(1)=subplot(211);hold on; box on
plot(BIS_accessible.time,BIS_accessible.signals.values,'.-','linewidth',2);
 yline(info.desired_BIS)
yline(info.upperLimit_BIS, 'r-'); yline(info.lowerLimit_BIS, 'r-')
xline(60, 'b--')
xlabel('Time [s]');  title('BIS - MPC Controller')
xlim([0 info.Tsim/6]); ylim([10 100]); set(gca, 'fontsize', 16); 

ax(2)=subplot(212); hold on; box on
plot(RASS_accessible.time,RASS_accessible.signals.values,'.-','linewidth',2);
 yline(info.desired_RASS)
yline(info.upperLimit_RASS, 'r-'); yline(info.lowerLimit_RASS, 'r-')
xline(60, 'b--')
xlabel('Time [s]'); title('RASS - MPC Controller')
xlim([0 info.Tsim/6]); ylim([-6 0]);  set(gca, 'fontsize', 16); 

linkaxes(ax, 'x')
% Figure: State comparison (Kalman vs Accessible)
figure
sgtitle('State comparison')
for i = 1:10
    subplot(5, 2, i)
    hold on; box on
plot(KalmanState_dx.time,KalmanState_dx.signals.values(:,i), 'r'); 
plot(accesibleState_dx.time,accesibleState_dx.signals.values(:,i), '--b');
xlabel('Time [s]');  title(['Comparison between States ',num2str(i)])
xlim([0 info.Tsim/6]); ylim([-50 10]); 
end
legend('Kalman estimate state', 'Accessible state')
     
% FIGURE: MPC controller w/ Kalman Filter
ax=[];
figure('color', 'w'); 
sgtitle('State estimated with Kalman filter')
ax(1)=subplot(211);hold on; box on
plot(BIS_kalman.time,BIS_kalman.signals.values, '.-', 'linewidth', 2); 
yline(info.desired_BIS)
yline(info.upperLimit_BIS, 'r-'); yline(info.lowerLimit_BIS, 'r-')
xline(60, 'b--')
xlabel('Time [s]');  title('BIS ')
xlim([0 info.Tsim]); ylim([10 100]); set(gca, 'fontsize', 16); 

ax(2)=subplot(212); hold on; box on
plot(RASS_kalman.time,RASS_kalman.signals.values, '.-', 'linewidth', 2); 
yline(info.desired_RASS)
yline(info.upperLimit_RASS, 'r-'); yline(info.lowerLimit_RASS, 'r-')
xline(60, 'b--')
xlabel('Time [s]'); title('RASS - MPC Controller')
xlim([0 info.Tsim/6]); ylim([-6 0]);  set(gca, 'fontsize', 16); 

linkaxes(ax, 'x')
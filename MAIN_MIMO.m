%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%-----------------------Simulator inspired to:--------------------------%%
%"An Open Source Patient Simulator for Design and Evaluation of Computer 
%Based Multiple Drug Dosing Control for Anesthetic and Hemodynamic Variables"
% C. Ionescu, M. Neckebroek, Mi. Ghita and D. Copot
% IEEE Access; 2021; doi 10.1109/ACCESS.2021.3049880  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                

%% Reset
clearvars; clear memory;  clc; warning off
close all;
addpath(genpath('src'));
load('continuousModel.mat')
load('operatingPoint.mat')

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
wantNoisyNBM = 1;               % 0= I do NOT want noise on NBM, 1= I want noise on NBM

info.T_induction=0;
%BIS
info.initial_BIS=100;
info.desired_BIS = 50;
info.upperLimit_BIS= 60;
info.lowerLimit_BIS= 40;

%RASS
info.initial_RASS=0;
info.desired_RASS = -4;
info.upperLimit_RASS= -3;
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

q=1*eye(2); 
% r=[1/25 0;0 1/6.25];
% q= 0.1*eye(2);
r= 1*eye(2);
PH=100;
F= -[1 0; 0 1];
% F = zeros(2);
f= u_eq;
% f=[0;0];
% 
condDiscMat=...
BuildCondendsedMPCmatrices(discreteMat,r,q,F,f,PH);


%% KALMAN FILTER
Q_KF = diag((0.1*x_eq + min(x_eq)).^2);
R_KF = 0.001*eye(2);

% %Setting the new parameter for the MPC
q_Kalman=10*eye(2); 
% q_Kalman = [50 0; 0 10];
r_Kalman= 0.001*eye(2);
% PH=200;
% F= -[1 0; 0 1];
% f= u_eq;
% 
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



%FIGURE: Comparison between estimated state and accessible state

% figure('color', 'w'); 
% hold on; box on
% plot(KalmanState_dx.time,KalmanState_dx.signals.values, 'r'); 
% plot(accesibleState_dx.time,accesibleState_dx.signals.values, '--b');
% xlabel('Time [s]');  title('Comparison between States (1)')
% xlim([0 info.Tsim/6]); ylim([-100 100]); set(gca, 'fontsize', 16);

figure
sgtitle('State comparison')
for i = 1:10
    subplot(5, 2, i)
    hold on; box on
plot(KalmanState_dx.time,KalmanState_dx.signals.values(:,i), 'r'); 
plot(accesibleState_dx.time,accesibleState_dx.signals.values(:,i), '--b');
xlabel('Time [s]');  title(['Comparison between States ',num2str(i)])
xlim([0 info.Tsim/6]); ylim([-2 2]); 
end
legend('Kalman estimate state', 'Accessible state')
     
% FIGURE: MPC controller - Kalman Filter
ax=[];
figure('color', 'w'); 
sgtitle('State estimated with Kalman filter')
ax(1)=subplot(211);hold on; box on
plot(BIS_kalman.time,BIS_kalman.signals.values, '.-', 'linewidth', 2); 
yline(info.desired_BIS)
yline(info.upperLimit_BIS, 'r-'); yline(info.lowerLimit_BIS, 'r-')
xline(60, 'b--')
xlabel('Time [s]');  title('BIS ')
xlim([0 info.Tsim/6]); ylim([10 100]); set(gca, 'fontsize', 16); 

ax(2)=subplot(212); hold on; box on
plot(RASS_kalman.time,RASS_kalman.signals.values, '.-', 'linewidth', 2); 
yline(info.desired_RASS)
yline(info.upperLimit_RASS, 'r-'); yline(info.lowerLimit_RASS, 'r-')
xline(60, 'b--')
xlabel('Time [s]'); title('RASS - MPC Controller')
xlim([0 info.Tsim/6]); ylim([-6 0]);  set(gca, 'fontsize', 16); 


linkaxes(ax, 'x')

% %% FIGURE TOT
% ax=[];
% figure('color', 'w'); 
% ax(1)=subplot(221);hold on; box on
% plot(BIS_output, '.-', 'linewidth', 2); yline(50)
% xlabel('Time [s]');  title('BIS')
% xlim([0 info.Tsim]); ylim([10 100]); set(gca, 'fontsize', 16); 
% 
% ax(2)=subplot(222); hold on; box on
% plot(RASS_output, '.-', 'linewidth', 2); yline(-4)
% xlabel('Time [s]'); title('RASS')
% xlim([0 info.Tsim]); ylim([-6 0]);  set(gca, 'fontsize', 16); 
% 
% 
% ax(6)=subplot(223); hold on; box on
% plot(u_Propofol.time,u_Propofol.signals.values,'.-', 'linewidth', 2)
% xlabel('Time [s]'); ylabel('mg/(kg min)'); title('Propofol')
% xlim([0 info.Tsim]);  set(gca, 'fontsize', 16); 
% 
% ax(7)=subplot(224); hold on; box on
% plot(u_Remifentanil.time,u_Remifentanil.signals.values,'.-', 'linewidth', 2)
% xlabel('Time [s]'); ylabel('mg/(kg min)'); title('Remifentanil')
% xlim([0 info.Tsim]);  set(gca, 'fontsize', 16); 
% 
% linkaxes(ax, 'x')
% 
% %% FIGURE BIS - LIN
% ax=[];
% figure('color', 'w'); 
% ax(1)=subplot(211);hold on; box on
% plot(u_Propofol_cont.time,u_Propofol_cont.signals.values,'.-', 'linewidth', 2)
% xlabel('Time [s]'); ylabel('mg/(kg min)'); title('Propofol - Continuous model')
% xlim([0 info.Tsim]);  set(gca, 'fontsize', 16); 
% 
% ax(2)=subplot(212);hold on; box on
% plot(BIS_cont.time,BIS_cont.signals.values, '.-', 'linewidth', 2); 
% yline(info.desired_BIS)
% yline(info.upperLimit_BIS, 'r-'); yline(info.lowerLimit_BIS, 'r-')
% xlabel('Time [s]'); ylabel('BIS'); title('BIS')
% xlim([0 info.Tsim]); ylim([10 100]); set(gca, 'fontsize', 16); 
% linkaxes(ax, 'x')
% 
% %% FIGURE RASS
% ax=[];
% figure('color', 'w'); 
% ax(1)=subplot(211);hold on; box on
% plot(u_Remifentanil_cont.time,u_Remifentanil_cont.signals.values,'.-', 'linewidth', 2)
% xlabel('Time [s]'); ylabel('mg/(kg min)'); title('Remifentanil - Continuous model')
% xlim([0 info.Tsim]);  set(gca, 'fontsize', 16); 
% 
% ax(2)=subplot(212);hold on; box on
% plot(RASS_cont.time,RASS_cont.signals.values, '.-', 'linewidth', 2); 
% yline(info.desired_RASS)
% yline(info.upperLimit_RASS, 'r-'); yline(info.lowerLimit_RASS, 'r-')
% xlabel('Time [s]'); title('RASS')
% xlim([0 info.Tsim]); ylim([-6 0]);  set(gca, 'fontsize', 16); 
% linkaxes(ax, 'x')
% 
% %% FIGURE TOT - LIN SYS
% ax=[];
% figure('color', 'w'); 
% ax(1)=subplot(221);hold on; box on
% plot(BIS_cont.time,BIS_cont.signals.values, '.-', 'linewidth', 2); yline(50)
% xlabel('Time [s]');  title('BIS Cont')
% xlim([0 info.Tsim]); ylim([10 100]); set(gca, 'fontsize', 16); 
% 
% ax(2)=subplot(222); hold on; box on
% plot(RASS_cont.time,RASS_cont.signals.values, '.-', 'linewidth', 2); yline(-4)
% xlabel('Time [s]'); title('RASS Cont')
% xlim([0 info.Tsim]); ylim([-6 0]);  set(gca, 'fontsize', 16); 
% 
% 
% ax(6)=subplot(223); hold on; box on
% plot(u_Propofol_cont.time,u_Propofol_cont.signals.values,'.-', 'linewidth', 2)
% xlabel('Time [s]'); ylabel('mg/(kg min)'); title('Propofol')
% xlim([0 info.Tsim]);  set(gca, 'fontsize', 16); 
% 
% ax(7)=subplot(224); hold on; box on
% plot(u_Remifentanil_cont.time,u_Remifentanil_cont.signals.values,'.-', 'linewidth', 2)
% xlabel('Time [s]'); ylabel('mg/(kg min)'); title('Remifentanil')
% xlim([0 info.Tsim]);  set(gca, 'fontsize', 16); 
% 
% linkaxes(ax, 'x')
% 
% 
% %% FIGURE COMPARISON
% 
% figure('color', 'w'); 
% ax(1)=subplot(211);hold on; box on
% plot(BIS_output, '.-', 'linewidth', 2); 
% plot(BIS_cont.time,BIS_cont.signals.values, '.-', 'linewidth', 2); yline(50)
% xlabel('Time [s]');  title('BIS')
% xlim([0 info.Tsim]); ylim([10 100]); set(gca, 'fontsize', 16); 
% legend('non-linear','linearized')
% 
% ax(2)=subplot(212); hold on; box on
% plot(RASS_output, '.-', 'linewidth', 2);
% plot(RASS_cont.time,RASS_cont.signals.values, '.-', 'linewidth', 2); yline(50)
% xlabel('Time [s]'); title('RASS')
% xlim([0 info.Tsim]); ylim([-6 0]);  set(gca, 'fontsize', 16); 
% legend('non-linear','linearized')
% 
% linkaxes(ax, 'x')
% 
% 
% 
% figure('color', 'w'); 
% ax(1)=subplot(211);hold on; box on
% plot(BIS_noNoise.time,BIS_noNoise.signals.values, '.-', 'linewidth', 2); 
% plot(BIS_cont.time,BIS_cont.signals.values, '.-', 'linewidth', 2); yline(50)
% xlabel('Time [s]');  title('BIS')
% xlim([0 info.Tsim]); ylim([10 100]); set(gca, 'fontsize', 16); 
% legend('non-linear, noise free','linearized')
% 
% ax(2)=subplot(212); hold on; box on
% plot(RASS_withoutDiscretization.time,RASS_withoutDiscretization.signals.values, '.-', 'linewidth', 2);
% plot(RASS_cont.time,RASS_cont.signals.values, '.-', 'linewidth', 2); yline(50)
% xlabel('Time [s]'); title('RASS')
% xlim([0 info.Tsim]); ylim([-6 0]);  set(gca, 'fontsize', 16); 
% legend('non-linear, non-discretized','linearized')
% 
% linkaxes(ax, 'x')
% 
% 
% 
% figure('color', 'w'); 
% ax(1)=subplot(211);hold on; box on
% plot(BIS_output, '.-', 'linewidth', 2); 
% plot(BIS_noNoise.time,BIS_noNoise.signals.values, '.-', 'linewidth', 2); yline(50)
% xlabel('Time [s]');  title('BIS')
% xlim([0 info.Tsim]); ylim([10 100]); set(gca, 'fontsize', 16); 
% legend('BIS with noise','BIS without noise')
% 
% ax(2)=subplot(212); hold on; box on
% plot(RASS_output, '.-', 'linewidth', 2);
% plot(RASS_withoutDiscretization.time,RASS_withoutDiscretization.signals.values, '.-', 'linewidth', 2); yline(50)
% xlabel('Time [s]'); title('RASS')
% xlim([0 info.Tsim]); ylim([-6 0]);  set(gca, 'fontsize', 16); 
% legend('RASS after discretization','RASS before discretization')
% 
% linkaxes(ax, 'x')
% 
% %% FIGURE TOT - LIN SYS DISC
% ax=[];
% figure('color', 'w'); 
% ax(1)=subplot(221);hold on; box on
% plot(BIS_linDisc.time,BIS_linDisc.signals.values, '.-', 'linewidth', 2); yline(50)
% xlabel('Time [s]');  title('BIS Disc')
% xlim([0 info.Tsim]); ylim([10 100]); set(gca, 'fontsize', 16); 
% 
% ax(2)=subplot(222); hold on; box on
% plot(RASS_linDisc.time,RASS_linDisc.signals.values, '.-', 'linewidth', 2); yline(-4)
% xlabel('Time [s]'); title('RASS Disc')
% xlim([0 info.Tsim]); ylim([-6 0]);  set(gca, 'fontsize', 16); 
% 
% 
% ax(6)=subplot(223); hold on; box on
% plot(u_Propofol_disc.time,u_Propofol_disc.signals.values,'.-', 'linewidth', 2)
% xlabel('Time [s]'); ylabel('mg/(kg min)'); title('Propofol')
% xlim([0 info.Tsim]);  set(gca, 'fontsize', 16); 
% 
% ax(7)=subplot(224); hold on; box on
% plot(u_Remifentanil_disc.time,u_Remifentanil_disc.signals.values,'.-', 'linewidth', 2)
% xlabel('Time [s]'); ylabel('mg/(kg min)'); title('Remifentanil')
% xlim([0 info.Tsim]);  set(gca, 'fontsize', 16); 
% 
% linkaxes(ax, 'x')
% 
% %% COMPARISON BETWEEN NON-LIN, LIN-CONT AND LIN-DISC
% 
% figure; 
% subplot(211);hold on; box on
% plot(BIS_output,'b'); 
% plot(BIS_cont.time,BIS_cont.signals.values,'r','LineWidth',2); 
% plot(BIS_linDisc.time, BIS_linDisc.signals.values,'g','LineWidth',2);
% yline(50)
% xlabel('Time [s]');  title('BIS')
% xlim([0 info.Tsim]); ylim([10 100]); set(gca, 'fontsize', 16); 
% legend('non-linear','linearized continuous','linearized discrete')
% 
% subplot(212);hold on; box on
% plot(RASS_output, 'b'); 
% plot(RASS_cont.time,RASS_cont.signals.values,'r', 'LineWidth', 2); 
% plot(RASS_linDisc.time, RASS_linDisc.signals.values,'g', 'LineWidth', 2);
% yline(50)
% xlabel('Time [s]');  title('RASS')
% xlim([0 info.Tsim]); ylim([-6 0]); set(gca, 'fontsize', 16); 
% legend('non-linear','linearized continuous','linearized discrete')
% 

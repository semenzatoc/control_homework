function [propofol, remifentanil, RASS, NMB] = drugs_parameters(patient)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%% Transfer function (PK) from Propofol input to effect site concentration   
V1p = 4.27; V2p = 18.9 - 0.391*(patient.age-53); V3p = 238; %[l]
Cl1p = 1.89 + 0.0456*(patient.weight - 77) - 0.0681*(patient.lbm - 59)+ 0.0264*(patient.height - 177);
Cl2p = 1.29 - 0.024*(patient.age - 53); Cl3p = 0.836; % [l.min^(-1)]
k10p = Cl1p/V1p; k12p = Cl2p/V1p; k13p = Cl3p/V1p; 
k21p = Cl2p/V2p; k31p = Cl3p/V3p; %[min^(-1)]
k1ep = 0.456; ke0p = 0.456; %[min^(-1)]

%% Hypnosis
%Input is Propofol, output is BIS (bispectral index)
%Transfer function from Propofol input to effect site concentration
Ap = [-(k10p + k12p + k13p) k21p k31p 0; k12p -k21p 0 0; k13p 0 -k31p 0; k1ep 0 0 -ke0p];
Bp = [1;0;0;0];
Cp = [0 0 0 1];
Dp = 0;
propofol.propSS = ss(Ap,Bp,Cp,Dp);

%% Constant declaration anesthetic models
hill_Propofol.gammaP = 2; hill_Propofol.C50P = 4.16; %Hill function Propofol
hill_Propofol.beta = 2.13; hill_Propofol.gamma = 4; hill_Propofol.sigma = 8.2; 
hill_Propofol.V1p = V1p;
propofol.hill_Propofol= hill_Propofol;

hill_Remifentanil.gammaR = 2.4; hill_Remifentanil.C50R = 8.84; %Hill function Remifentanil
remifentanil.hill_Remifentanil = hill_Remifentanil;

%% Transfer function (PK) from Remifentanil input to effect site concentration
V1r = 5.1 - 0.0201*(patient.age-40) + 0.072*(patient.lbm-55); 
V2r = 9.82 - 0.0811*(patient.age-40) + 0.108*(patient.lbm-55); V3r = 5.42; %[l]
Cl1r = 2.6 + 0.0162*(patient.age - 40) + 0.0191*(patient.lbm-55); %[l*min^(-1)]
Cl2r = 2.05 - 0.0301*(patient.age - 40); Cl3r = 0.076 - 0.00113*(patient.age - 40); %[l*min^(-1)]
k10r = Cl1r/V1r; k12r = Cl2r/V1r; k13r = Cl3r/V1r; k21r = Cl2r/V2r; k31r = Cl3r/V3r; %[min^(-1)]
ke0r = 0.595 - 0.007*(patient.age - 40); k1er = 0.456; %[min^(-1)]
remifentanil.hill_Remifentanil.V1r = V1r;

%% Transfer function from Remifentanil effect site concentration to level RASS sedation score
remifentanil.k1r = 0.81; remifentanil.k0r = 0.81;

%% Transfer function from Atracurium input to effect site concentration
k1a = 1; k2a = 4; k3a = 10; 
NMB.alphaNMB = 0.0374; NMB.gammaN = 2.6677; NMB.C50N = 3.2425;

%% %% Interaction Prop to hemodynamic
%%%============================================
%Interaction model from Propofol to CO
propofol.k1Pco = 0.81; propofol.k0Pco = 0.81; propofol.E0Pco = 5; 
propofol.EmaxPco = 5; propofol.gainPco=10; propofol.gammaPco = 4.5; propofol.C50Pco = 8;

%Interaction model from Propofol to MAP
propofol.k1Pmap = 0.61; propofol.k0Pmap = 0.81; propofol.E0Pmap = 5; 
propofol.EmaxPmap = 5; propofol.gainPmap=15; propofol.gammaPmap = 4.5; propofol.C50Pmap = 6;


%% Interaction Remi to hemodynamic
%Interaction model from Remi to CO
remifentanil.k1Rco = 0.51; remifentanil.k0Rco = 0.51; remifentanil.E0Rco = 15; 
remifentanil.EmaxRco = 5; remifentanil.gainRco=10; remifentanil.gammaRco = 4.5; 
remifentanil.C50Rco = 12;

%Interaction model from Remi to MAP
remifentanil.k1Rmap = 0.31; remifentanil.k0Rmap = 0.31; remifentanil.E0Rmap = 70; 
remifentanil.EmaxRmap = 70; remifentanil.gainRmap=10; remifentanil.gammaRmap = 4.5; 
remifentanil.C50Rmap = 17;





%% Analgesia
%Transfer function from Remifentanil input to effect concentration
Ar = [-(k10r + k12r + k13r) k21r k31r 0; k12r -k21r 0 0; k13r 0 -k31r 0; k1er 0 0 -ke0r];
Br = [1;0;0;0];
Cr = [0 0 0 1];
Dr = 0;
remiSS = ss(Ar,Br,Cr,Dr);
remifentanil.remiSS = remiSS;

%% Model from Remifentanil to level RASS sedation score
sysRASS = tf(-2, [1 2]); %Generate transfer function of continuous model
[RASS.numRass, RASS.denRass]=tfdata(sysRASS,'v'); %Get coefficients of transfer function of continuous model
% Calculate the gain for RASS
RASS.dcgRASS = dcgain(tf(RASS.numRass,RASS.denRass)*tf(1,[remifentanil.k1r*15 remifentanil.k0r]));

%% Neuromuscular blockade (NMB)
% NMB model
sysNMB=zpk([],[-k1a*NMB.alphaNMB, -k2a*NMB.alphaNMB, -k3a*NMB.alphaNMB],...
    k1a*k2a*k3a*NMB.alphaNMB^3); %Generate transfer function of nominal model
[numN,denN]=tfdata(sysNMB,'v'); %Get coefficients of transfer function of nominal model
[Anmb,Bnmb,Cnmb,Dnmb] = tf2ss(numN,denN); %Convert transfer function to state space
NMB.nmbSS = ss(Anmb,Bnmb,Cnmb,Dnmb); %Generate state space model




end


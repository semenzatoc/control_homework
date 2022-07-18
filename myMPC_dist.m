function u_out = myMPC_dist(input,condDiscMat,PH)
ref_bis = input(1); 
ref_rass = input(2);
input_prop = input(3);
input_remi = input(4);
% input vector
u_call = repmat([input_prop; input_remi],PH,1)';
% state vector
x = input(5:14);
% y reference vector
% ref = repmat([50; -4],PH,1);
ref = repmat([ref_bis; ref_rass],PH,1);
%disturbance
d = zeros(2*PH,1);

AC =condDiscMat.call_AC;
BC =condDiscMat.call_BC;
MC =condDiscMat.call_MC;
QC =condDiscMat.call_Q;
RC =condDiscMat.call_R;
FC =condDiscMat.call_F;
fC =condDiscMat.call_f;

Q = BC'*QC*BC + RC;
c = ((AC*x + MC*d - ref)'*QC*BC- u_call*RC)';
A = FC;
b = fC;

options = optimset('Display','off','Algorithm','interior-point-convex');
% u_out= quadprog(Q,c,A,b);
[u_sequence,~,~]=quadprog((Q+Q')/2,c,A,b,[],[],[],[],[],options);

u_out(1)= u_sequence(1);
u_out(2)= u_sequence(2);



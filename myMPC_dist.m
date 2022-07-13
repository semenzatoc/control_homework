function u_out = myMPC_dist(input,condDiscMat,PH)
input_prop = input(1);
input_remi = input(2);
% input vector
u_call =repmat([input_prop; input_remi],PH,1)';
% state vector
x = input(3:12);
% y reference vector
% ref = repmat([50; -4],PH,1);
ref = zeros(2*PH,1);
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



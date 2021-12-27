load moffed_cropped.mat

[P,M]=size(Y)
Px=50;
r=3;
q=0.5;
alpha=0.5;
maxG=1;
N=r*(r-1)/2;



ascZ=0;
delta=1;
% Initial values for endmenbers
Avca=vca(Y','endmembers',r);
Svca= hyperFcls(Y',Avca)';
Bvca=createB(Avca);

Bd=createB(Avca);
Bd=[Bd; zeros(1,N)];

Ad=[Avca;delta*ones(1,r)];
Yd=[Y delta*ones(P,1)];


% Initial values for abundances
Z=nonegReg((Yd-Svca*Ad')',Bd,[])';
Z=0.1*rand(size(Z));
hZ=0.001;

% Estimate parameters
[hS, all_hss, all_ebic]=estimate_bil_hs_WithEbic(Yd,Ad,Bd,Svca,Z,q,alpha,hZ);
[hZ, all_hzs, all_ebic]=estimate_bil_hz_WithEbic(Yd,Ad,Bd,Svca,Z,q,alpha,hS);

% Run bilinear method
[Abil,Bbil,Sbil,Zbil]=sparseBilinearUnmixing(Y,Avca,Svca,Z,hS,hZ,delta,q,ascZ);


% Fan method
Sb0=createB(Svca);
Bd=createB(Avca);
Bd=[Bd; zeros(1,N)];

G0=maxG*rand(P,N);
G0(Sb0==0)=0;
[lambda_fan, all_lambdas, all_ebic]=estimateLambdaWithEbicFan(Yd,Ad,Bd,Svca,maxG,q,alpha,0.02);

[A,B,S,G]=sparseBilinearFanUnmixing(Y,Avca,Svca,G0,maxG,lambda_fan,delta,q);
Z=createB(S);









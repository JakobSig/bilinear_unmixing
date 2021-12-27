function [lambda, all_lambdas, all_ebic]=estimateLambdaWithEbicFan(Y,Ad,Bd,S,maxG,q,alpha,d)
%(Y,Ad,Bd,S,hS,hZ,q,stopCriteria,checkConv)
% Estimate lambda with EBIC
done=0;
%lambdas=[0.001 0.005 0.01 0.015];
%lambdas=[0.0010    0.0023    0.0037    0.0050]*3;

Yo=Y(:,1:end-1);
Ao=Ad(1:end-1,:);


%d=0.02;
lambdas=[0.0010 d 2*d 3*d];

ebic_iter=10000;
checkConv=10;
stopCriteria=10^-4;
[M,r]=size(Ad);
%A0=[Ad Bd];
delta=10;
Z=createB(S);
%u=10*norm(Y)/(numel(Y));
u=0.01*norm(Y,'fro').^2/(numel(Y));
G=sparseBilinearUnmixingFanGStepSparsePen(Yo,Ao,maxG*rand(size(Z)),S,Z,stopCriteria,checkConv,maxG,u,q,0,delta);

for i=1:4
    S0{i}=sparseBilinearUnmixingFanSStep(Y,Ad,Bd,G,S,lambdas(i),q,delta,stopCriteria,checkConv);
    Z=createB(S0{i});
    G0=sparseBilinearUnmixingFanGStepSparsePen(Yo,Ao,G,S0{i},Z,stopCriteria,checkConv,maxG,u,q,lambdas(i),delta);
    ebic(i)=EBICfan(Yo,S0{i},Ao,G0,alpha);
end
[a,idx]=min(ebic);
iter=0;
all_lambdas=[lambdas];all_ebic=[ebic];
while done==0
    iter=iter+1;
    if idx==1
        lambdas(1)=max(lambdas(1)-d,0);
        lambdas(4)=lambdas(2); ebic(4)=ebic(2);
        d=(lambdas(4)-lambdas(1))/3;
        lambdas(2)=lambdas(1)+d;
        lambdas(3)=lambdas(1)+2*d;
        for i=1:3
            S0{i}=sparseBilinearUnmixingFanSStep(Y,Ad,Bd,G,S,lambdas(i),q,delta,stopCriteria,checkConv);
            Z=createB(S0{i});
            G0=sparseBilinearUnmixingFanGStepSparsePen(Yo,Ao,G,S0{i},Z,stopCriteria,checkConv,maxG,u,q,lambdas(i),delta);
            ebic(i)=EBICfan(Yo,S0{i},Ao,G0,alpha);
        end
    end
    if idx==2
        lambdas(4)=lambdas(3); ebic(4)=ebic(3);
        d=(lambdas(4)-lambdas(1))/3;
        lambdas(2)=lambdas(1)+d;
        lambdas(3)=lambdas(1)+2*d;
        for i=[2 3]
            S0{i}=sparseBilinearUnmixingFanSStep(Y,Ad,Bd,G,S,lambdas(i),q,delta,stopCriteria,checkConv);
            Z=createB(S0{i});
            G0=sparseBilinearUnmixingFanGStepSparsePen(Yo,Ao,G,S0{i},Z,stopCriteria,checkConv,maxG,u,q,lambdas(i),delta);
            ebic(i)=EBICfan(Yo,S0{i},Ao,G0,alpha);
        end
    end
    if idx==3
        lambdas(1)=lambdas(2);ebic(1)=ebic(2);S0{1}=S0{2};
        d=(lambdas(4)-lambdas(1))/3;
        lambdas(2)=lambdas(1)+d;
        lambdas(3)=lambdas(1)+2*d;
        for i=[2 3]
            S0{i}=sparseBilinearUnmixingFanSStep(Y,Ad,Bd,G,S,lambdas(i),q,delta,stopCriteria,checkConv);
            Z=createB(S0{i});
            G0=sparseBilinearUnmixingFanGStepSparsePen(Yo,Ao,G,S0{i},Z,stopCriteria,checkConv,maxG,u,q,lambdas(i),delta);
            ebic(i)=EBICfan(Yo,S0{i},Ao,G0,alpha);
        end
    end
    if idx==4
        lambdas(1)=lambdas(3);ebic(1)=ebic(3);S0{1}=S0{3};
        lambdas(2)=lambdas(4);ebic(2)=ebic(4);
        lambdas(3)=lambdas(1)+2*d;
        lambdas(4)=lambdas(1)+3*d;
        for i=3:4
            S0{i}=sparseBilinearUnmixingFanSStep(Y,Ad,Bd,G,S,lambdas(i),q,delta,stopCriteria,checkConv);
            Z=createB(S0{i});
            G0=sparseBilinearUnmixingFanGStepSparsePen(Yo,Ao,G,S0{i},Z,stopCriteria,checkConv,maxG,u,q,lambdas(i),delta);
            ebic(i)=EBICfan(Yo,S0{i},Ao,G0,alpha);
        end
    
        
    end
    [a,idx]=min(ebic);
    
    all_lambdas=[all_lambdas lambdas];
    all_ebic=[all_ebic ebic];
    [a,b]=sort(all_lambdas);
    figure(1);hold off;plot(all_lambdas(b),all_ebic(b),'x-'); hold on;plot(lambdas,ebic,'rx-');hold off;shg;drawnow
    if abs((lambdas(1)-lambdas(4))/lambdas(4)) < 0.1 || iter>50
        done=1
        lambda=lambdas(idx);
    end
    %A=A0{1};S=S0{1};
end
[a,b]=min(all_ebic);
lambda=all_lambdas(b)

end
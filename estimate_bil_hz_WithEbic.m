function [hz, all_hzs, all_ebic]=estimate_bil_hz_WithEbic(Y,Ad,Bd,S,Z,q,alpha,hs)
%(Y,Ad,Bd,S,hS,hZ,q,stopCriteria,checkConv)
% Estimate hs with EBIC
done=0;
%hzs=[0.001 0.005 0.01 0.015];
S=[S Z];
[M,r]=size(Ad);
[M,N]=size(Bd);
A0=[Ad Bd];

scaler=5;
hzs=[0.0010    0.0023    0.0037    0.0050]*scaler;
hs=repmat(hs,1,r);
d=0.002*scaler;
ebic_iter=10000;
checkConv=10;
stopCriteria=10^-5;


for i=1:4
    S0{i}=sparseBilinearUnmixingSStep(Y,Ad,Bd,S,hs,repmat(hzs(i),1,N),q,stopCriteria,checkConv);
    ebic(i)=EBIC(Y,S0{i},A0,alpha);
end
[a,idx]=min(ebic);
iter=0;
all_hzs=[hzs];all_ebic=[ebic];

while done==0
    iter=iter+1;
    if idx==1
        hzs(1)=max(hzs(1)-d,0);
        hzs(4)=hzs(2); ebic(4)=ebic(2);
        d=(hzs(4)-hzs(1))/3;
        hzs(2)=hzs(1)+d;
        hzs(3)=hzs(1)+2*d;
        for i=1:3
            S0{i}=sparseBilinearUnmixingSStep(Y,Ad,Bd,S,hs,repmat(hzs(i),1,N),q,stopCriteria,checkConv);
            ebic(i)=EBIC(Y,S0{i},A0,alpha);
        end
    end
    if idx==2
        hzs(4)=hzs(3); ebic(4)=ebic(3);
        d=(hzs(4)-hzs(1))/3;
        hzs(2)=hzs(1)+d;
        hzs(3)=hzs(1)+2*d;
        for i=[2 3]
            S0{i}=sparseBilinearUnmixingSStep(Y,Ad,Bd,S,hs,repmat(hzs(i),1,N),q,stopCriteria,checkConv);
            ebic(i)=EBIC(Y,S0{i},A0,alpha);
        end
    end
    if idx==3
        hzs(1)=hzs(2);ebic(1)=ebic(2);S0{1}=S0{2};
        d=(hzs(4)-hzs(1))/3;
        hzs(2)=hzs(1)+d;
        hzs(3)=hzs(1)+2*d;
        for i=[2 3]
            S0{i}=sparseBilinearUnmixingSStep(Y,Ad,Bd,S,hs,repmat(hzs(i),1,N),q,stopCriteria,checkConv);
            ebic(i)=EBIC(Y,S0{i},A0,alpha);
        end
    end
    if idx==4
        hzs(1)=hzs(3);ebic(1)=ebic(3);S0{1}=S0{3};
        hzs(2)=hzs(4);ebic(2)=ebic(4);
        hzs(3)=hzs(1)+2*d;
        hzs(4)=hzs(1)+3*d;
        for i=3:4
            S0{i}=sparseBilinearUnmixingSStep(Y,Ad,Bd,S,hs,repmat(hzs(i),1,N),q,stopCriteria,checkConv);
            ebic(i)=EBIC(Y,S0{i},A0,alpha);
        end
    
        
    end
    [a,idx]=min(ebic);
    
    all_hzs=[all_hzs hzs];
    all_ebic=[all_ebic ebic];
    [a,b]=sort(all_hzs);
    figure(1);plot(all_hzs(b),all_ebic(b),'x-'); hold on;plot(hzs,ebic,'rx-');hold off;shg;drawnow
    if abs(hzs(1)-hzs(4))<0.004 || iter>20
        done=1
        hs=hzs(idx);
    end
    %A=A0{1};S=S0{1};
end
[a,b]=min(all_ebic);
hz=all_hzs(b)

end
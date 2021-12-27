function [hs, all_hss, all_ebic]=estimate_bil_hs_WithEbic(Y,Ad,Bd,S,Z,q,alpha,hZ)
%(Y,Ad,Bd,S,hS,hZ,q,stopCriteria,checkConv)
% Estimate hs with EBIC
done=0;
%hss=[0.001 0.005 0.01 0.015];
S=[S Z];
[M,r]=size(Ad);
[M,N]=size(Bd);
A0=[Ad Bd];

scaler=5;
hss=[0.0010    0.0023    0.0037    0.0050]*scaler;
hz=repmat(hZ,1,N);
d=0.002*scaler;
ebic_iter=10000;
checkConv=10;
stopCriteria=10^-5;


for i=1:4
    S0{i}=sparseBilinearUnmixingSStep(Y,Ad,Bd,S,repmat(hss(i),1,r),hz,q,stopCriteria,checkConv);
    ebic(i)=calcEbic(Y,S0{i},A0,r,alpha);%EBIC(Y,S0{i},A0,alpha);
end
[a,idx]=min(ebic);
iter=0;
all_hss=[hss];all_ebic=[ebic];

while done==0
    iter=iter+1;
    if idx==1
        hss(1)=max(hss(1)-d,0);
        hss(4)=hss(2); ebic(4)=ebic(2);
        d=(hss(4)-hss(1))/3;
        hss(2)=hss(1)+d;
        hss(3)=hss(1)+2*d;
        for i=1:3
            S0{i}=sparseBilinearUnmixingSStep(Y,Ad,Bd,S,repmat(hss(i),1,r),hz,q,stopCriteria,checkConv);
            ebic(i)=calcEbic(Y,S0{i},A0,r,alpha);%EBIC(Y,S0{i},A0,alpha);
        end
    end
    if idx==2
        hss(4)=hss(3); ebic(4)=ebic(3);
        d=(hss(4)-hss(1))/3;
        hss(2)=hss(1)+d;
        hss(3)=hss(1)+2*d;
        for i=[2 3]
            S0{i}=sparseBilinearUnmixingSStep(Y,Ad,Bd,S,repmat(hss(i),1,r),hz,q,stopCriteria,checkConv);
            ebic(i)=calcEbic(Y,S0{i},A0,r,alpha);%EBIC(Y,S0{i},A0,alpha);
        end
    end
    if idx==3
        hss(1)=hss(2);ebic(1)=ebic(2);S0{1}=S0{2};
        d=(hss(4)-hss(1))/3;
        hss(2)=hss(1)+d;
        hss(3)=hss(1)+2*d;
        for i=[2 3]
            S0{i}=sparseBilinearUnmixingSStep(Y,Ad,Bd,S,repmat(hss(i),1,r),hz,q,stopCriteria,checkConv);
            ebic(i)=calcEbic(Y,S0{i},A0,r,alpha);%EBIC(Y,S0{i},A0,alpha);
        end
    end
    if idx==4
        hss(1)=hss(3);ebic(1)=ebic(3);S0{1}=S0{3};
        hss(2)=hss(4);ebic(2)=ebic(4);
        hss(3)=hss(1)+2*d;
        hss(4)=hss(1)+3*d;
        for i=3:4
            S0{i}=sparseBilinearUnmixingSStep(Y,Ad,Bd,S,repmat(hss(i),1,r),hz,q,stopCriteria,checkConv);
            ebic(i)=calcEbic(Y,S0{i},A0,r,alpha);%EBIC(Y,S0{i},A0,alpha);
        end
    
        
    end
    [a,idx]=min(ebic);
    
    all_hss=[all_hss hss];
    all_ebic=[all_ebic ebic];
    [a,b]=sort(all_hss);
    figure(1);plot(all_hss(b),all_ebic(b),'x-'); hold on;plot(hss,ebic,'rx-');hold off;shg;drawnow
    if abs(hss(1)-hss(4))<0.004 || iter>20
        done=1
        hs=hss(idx);
    end
    %A=A0{1};S=S0{1};
end
[a,b]=min(all_ebic);
hs=all_hss(b)

end


function ebic=calcEbic(Y,S,A,r,alpha)
    ebic=EBIC(Y,S,A,alpha);
    %Z=S(:,r+1:end);
    %B=A(:,r+1:end);
    %S=S(:,1:r);
    %A=A(:,1:r);
    %ebic=EBIC(Y-Z*B',S,A,alpha);
end

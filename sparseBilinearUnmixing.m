function [A,B,S,Z]=sparseBilinearUnmixing(Y,A,S,Z,hS,hZ,delta,q,ascZ)
% Model: Y=SA^T+ZB^T
% Method solves: argmin_{S,A,Z,B} 1/2||Y-SA^T-ZB^T||^2_F + h_S ||S||_q + h_Z ||Z||_q
% Y data
% A endmembers
% S linear abundances
% B bilinear combos of endmembers
% Z bilinear abundances
% ascZ=1 if Z is to be included in the ASC, 0 otherwize 

    tic
    [M,r]=size(A);
    [P,~]=size(Y);
    maxIters=10000;
    stopCriteria=10^-5;
    mainStopCriteria=10^-4;

    B=createB(A);
    Nq=createNq(r);
    Sold=S;
    Aold=A;
    Zold=Z;
    N=r*(r-1)/2;
    hS=repmat(hS,1,r);
    hZ=repmat(hZ,1,N);
    
    deltaZ=0;
    if ascZ==1
      deltaZ=delta;  
    end
    
	Ad=[A;delta*ones(1,r)];
    Bd=[B;deltaZ*ones(1,N)];
    Yd=[Y delta*ones(P,1)];
    global checkConv
    for iter=1:maxIters
        
        checkConv=5;
        %disp(['Iteration ' num2str(iter)])
        fprintf('.')
        SZ=[S Z];
        SZ=sparseBilinearUnmixingSStep(Yd,Ad,Bd,SZ,hS,hZ,q,stopCriteria,checkConv);
        S=SZ(:,1:r);
        Z=SZ(:,r+1:end);
        
        % min_E ||Y-AS-MB||^2_F
        [A,B]=Astep_dyadic(Y,B,Z,S,A,stopCriteria,Nq,0);
        Ad=[A;delta*ones(1,r)];
        Bd=[B;deltaZ*ones(1,N)];
        if rem(iter,10)==0
            normA=norm(A-Aold,'fro')/norm(A,'fro');
            normS=norm(S-Sold,'fro')/norm(S,'fro');
            normZ=0.1*norm(Z-Zold,'fro')/max(10^-6,norm(Z,'fro'));
            %[normE normA normB]
            if(normA<mainStopCriteria && normS<mainStopCriteria && normZ<mainStopCriteria )
                break;
            end
        end
        Aold=A;Sold=S;Zold=Z;
        
    end
    toc
    
end






 function [A,B]=Astep_dyadic(Y,B,Z,S,A,stopCriteria,Nq,hA)
     global checkConv;

    [M,r]=size(A);
    indx=zeros(r,r-1);
    for i=1:r
        tmp=1:r;
        tmp(i)=[];
        indx(i,:)=tmp;
    end
    N=r*(r-1)/2;
    aN=1:N;
    for i=1:10000
        Aold=A;
        for q=1:r
            nq=Nq(q,:);
            neq=aN;neq(nq)=[];
            mj=indx(q,:);
            Xq=S(:,mj)*A(:,mj)'+Z(:,neq)*B(:,neq)';
            Yq=Y-Xq;
           
            Sq=S(:,q)*ones(M,1)';
            ZA=Z(:,nq)*A(:,mj)';
            Zq=Sq+ZA;
            aq=calcSlowEq(Yq,Zq,hA);
            aq=max(0,aq);
            A(:,q)=aq; 
            B=createB(A);
        end
         if rem(i,checkConv)==0
            normA=norm(A-Aold,'fro')/norm(A,'fro');
            if(normA<stopCriteria)
                return;
            end
            
        end
    end
    B=createB(A);
 end



 
 
 
function aq=calcSlowEq(YmX,Z,hA)
    [~,M]=size(YmX);
    aq=zeros(M,1);
    for m=1:M
        z=Z(:,m);
        aq(m)=YmX(:,m)'*z/(z'*z+hA);
    end
end






% Row q hold indexes in M that have either eq*ei or ei*eq
function Nq=createNq(r)
    N=r*(r-1)/2;
    Nq = tril(ones(r-1));
    Nq(Nq==1) = 1:N;
    Lq=zeros(r);
    Lq(2:end,1:end-1)=Nq;
    Uq=Lq';
    Nq = Lq(:, 1:end-1) + Uq(:, 2:end);
end



















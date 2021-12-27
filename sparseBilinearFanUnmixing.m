function [A,B,S,G]=sparseBilinearFanUnmixing(Y,A,S,G,maxG,hS,delta,q)
% Model: Y=SA^T+ZB^T
% Method solves: argmin_{S,A,Z,B} 1/2||Y-SA^T-ZB^T||^2_F + h_S ||S||_{q,1}
% s.t. A>0, S>0, 0<G<maxG,
% Y data (PxM) P=#pixels, M=#spectral bands
% A endmembers
% S linear abundances
% G Gamma - 0 < g < maxG
% B bilinear endmembers =[(A_1 * A_2) (A_1 * A_3) ... ]
% Z bilinear abundances Z=[ (G_1 * S_1 * S_2) (G_2 * S_1 * S_3) ... ]
    %Px=50;

    tic
    [M,r]=size(A);
    [P,N]=size(G);
    maxIters=10000;
    stopCriteria=10^3;%10^-5;
    mainStopCriteria=10^-4;

    if maxG==0
        G=zeros(P,N);
    end
    
    B=createB(A);
    Z=createB(S);
    Nq=createNq(r);
    Aold=A;Sold=S;Gold=G;
    %Z=createZ(G,S);
    N=r*(r-1)/2;
    %hS=repmat(hS,1,r);
    
    u=0.01*norm(Y,'fro').^2/(numel(Y));
    
	Ad=[A;delta*ones(1,r)];
    Bd=[B;zeros*ones(1,N)];
    Yd=[Y delta*ones(P,1)];
    global checkConv
    costIter=1;
    cols=[ 'm' 'c' 'r' 'b' 'k']; str=num2str(cputime); s =str2num(str(end))+1;c=ceil(s/2); plotcol=cols(c);
    for iter=1:maxIters
        checkConv=5;
        if maxG>0
            G=sparseBilinearUnmixingFanGStepSparsePen(Y,A,G,S,Z,stopCriteria,checkConv,maxG,u,q,hS,delta);
        end
        [S,Z]=sparseBilinearUnmixingFanSStep(Yd,Ad,Bd,G,S,hS,q,delta,stopCriteria,checkConv);
        [A,B]=Astep_dyadic(Y,B,G,S,A,stopCriteria/10,Nq,0,checkConv);
        Ad=[A;delta*ones(1,r)];
        Bd=[B;zeros*ones(1,N)];
        if rem(iter,10)==0
            normA=norm(A-Aold,'fro')/norm(A,'fro');
            normS=norm(S-Sold,'fro')/norm(S,'fro');
            normG=0.01*norm(G-Gold,'fro')/max(10^-6,norm(G,'fro'));
            %[normE normA normB]
            if(normA<mainStopCriteria && normS<mainStopCriteria && normG<mainStopCriteria )
                break;
            end
        end
        Aold=A;Sold=S;Gold=G;
        
    end
    toc
    
end

function [A,B]=Astep_mm(Y,B,S,G,Z,A,stopCriteria,checkConv)
    GZ=G.*Z;
    for i=1:10000
        Aold=A;
        X=max(0,Y-GZ*B');
        A=(A.*(X'*S))./(A*(S'*S));
        B=createB(A);
        if rem(i,checkConv)==0
            normA=norm(A-Aold,'fro')/norm(A,'fro');
            if(normA<stopCriteria)
                return;
            end
        end
    end
end


 function [A,B]=Astep_dyadic(Y,B,G,S,A,stopCriteria,Nq,hA,checkConv)
    Z=G.*createB(S);    
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
        for q=randperm(r)
            nq=Nq(q,:);
            neq=aN;neq(nq)=[];
            mj=indx(q,:);
            Xq=[S(:,mj) Z(:,neq)]*[A(:,mj) B(:,neq)]';
            Yq=Y-Xq;
           
            Zq=[S(:,q) Z(:,nq)]*[ones(M,1) A(:,mj)]';
            aq=calcAq(Yq,Zq);
            %aq=diag(Yq'*Zq)./sum(Zq.^2)';
            aq=max(0,aq);
            A(:,q)=aq; 
            B=updateB(A,B,q);
        end
         if rem(i,checkConv)==0
            normA=norm(A-Aold,'fro')/norm(A,'fro');
            if(normA<stopCriteria)
                return;
            end
            
        end
    end
 end



function aq=calcAq(YmX,Z)
    [~,M]=size(YmX);
    aq=zeros(M,1);
    for m=1:M
        z=Z(:,m);
        aq(m)=YmX(:,m)'*z/(z'*z);
    end
end
 
 
 
function aq=calcSlowEqhA(YmX,Z,hA)
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



















function G=sparseBilinearUnmixingFanGStepSparsePen(Y,A,G,S,Z,stopCriteria,checkConv,maxG,u,q,hS,delta)
    %stopCriteria=10^-5;
    [M,r]=size(A);
    G(Z==0)=0;
    [P,M]=size(Y);
    %nonZeros=logical(ones(size(Z)));
    nonZeros=(Z>0);
    B=createB(A);
    N=r*(r-1)/2;
    indx=zeros(N,N-1);
    %u=0.01;
  	for n=1:N
        tmp=1:N;
        tmp(n)=[];
        indx(n,:)=tmp;
    end
    %SAT=S*A';
    for n=1:N
        nonZ=nonZeros(:,n);
        Pk{n}=Z(nonZ,n)*B(:,n)';
        %[a,b]=hist(Pk{n}(:),100);
        BTB=B(:,n)'*B(:,n);
        tmp=Z(nonZ,n).^2*BTB+u;
        PkPk{n}=1./tmp;
    end
    
    YmSAT=Y-S*A';
    costIter=1;
    for Giter=1:10000
        Gold=G;
        for n=1:N
            nonZ=nonZeros(:,n);
            mn=indx(n,:);
            % Create Rk
            Rk=(G(nonZ,mn).*Z(nonZ,mn))*B(:,mn)';
            % Create Xk
            %Xk=Y(nonZ,:)-SAT(nonZ,:)-Rk;
            Xk=YmSAT(nonZ,:)-Rk;
            g=sum(Pk{n}.*Xk,2).*PkPk{n};
            g=max(0,g);
            G(nonZ,n)=min(maxG,g);
            %cost(costIter)=calcBilFanCost(Y,S,A,G,q,hS,delta,u);costIter=costIter+1;
        end
        %if max(diff(cost))>0 || any(isnan(cost))
        %    foo=0;
        %end
        if rem(Giter,checkConv)==0
            ndG=norm(G-Gold,'fro')/norm(Gold,'fro');
            if(ndG<stopCriteria)
                return;
            end
        end
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



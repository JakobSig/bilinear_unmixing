function [S,Z]=sparseBilinearUnmixingFanSStep(Y,A,B,G,S,hS,q,delta,stopCriteria,checkConv)

    [P,r]=size(S);
    %B=createB(A);
    N=r*(r-1)/2;
    Nk=createNq(r);
    aN=1:N;
  	for i=1:r
        tmp=1:r;
        tmp(i)=[];
        indx(i,:)=tmp;
    end
    nonZeros=S>0;
    pS=sum(nonZeros);
    
    hSq=(hS*q);
    qmI=q-1;
    %costIter=1;
    %cost(costIter)=calcBilFanCost(Y(:,1:end-1),S,A(1:end-1,:),G,q,hS,delta);costIter=costIter+1;
    Z=createB(S);
    %SAT=S*A';
    %multSum=@(a,b) sum(a.*b,2);
    for Siter=1:10000
        Sold=S;
        for k=1:r
            nonZ=nonZeros(:,k);
            nk=Nk(k,:);
            nek=aN;nek(nk)=[];
            mk=indx(k,:);
            
            Sn=S(nonZ,:);
            Gn=G(nonZ,:);

            
            % Create Rk
            %Rk=[S(:,mk) G(:,nek).*Z(:,nek)]*[A(:,mk) B(:,nek)]';
            Rk=[Sn(:,mk) Gn(:,nek).*Z(nonZ,nek)]*[A(:,mk) B(:,nek)]';
            % Create Pk
            %Pk=[G(:,nk).*S(:,mk) ones(P,1)]*[B(:,nk) A(:,k)]';
            Pk=[Gn(:,nk).*Sn(:,mk) ones(pS(k),1)]*[B(:,nk) A(:,k)]';

            %cost(costIter)=calcBilFanCost(Y(:,1:end-1),S,A(1:end-1,:),G,q,hS,delta);costIter=costIter+1;
            s=S(nonZ,k);
            if hS>0
                s=s.*sum(Pk.*(Y(nonZ,:)-Rk),2)./( s.*sum(Pk.*Pk,2)+hSq*exp(qmI*log(s)));
                %s=s.*sum(Pk.*(Y(nonZ,:)-Rk),2)./( s.*sum(Pk.*Pk,2)+hSq*s.^qmI);
            else
                %s=sum(Pk.*(Y-Rk),2)./sum(Pk.*Pk,2);
                s=sum(Pk.*(Y(nonZ,:)-Rk),2)./sum(Pk.*Pk,2);
            end
            %S(:,k)=max(0,s);
            S(nonZ,k)=max(0,s);
            Z=updateB(S,Z,k);
            %cost(costIter)=calcBilFanCost(Y(:,1:end-1),S,A(1:end-1,:),G,q,hS,delta);costIter=costIter+1;

            %for p=1:P
            %    s=S(p,k);
            %    s2(p)=(s*Pk(p,:)*X(p,:)')/(s*Pk(p,:)*Pk(p,:)'+(hS/q)*s^(q-1));
            %end
        end
        if rem(Siter,checkConv)==0
            ndS=norm(S-Sold,'fro')/norm(Sold,'fro');
            if(ndS<stopCriteria/100)
                return;
            end
        end
    end
end

function Nq=createNq(r)
    N=r*(r-1)/2;
    Nq = tril(ones(r-1));
    Nq(Nq==1) = 1:N;
    Lq=zeros(r);
    Lq(2:end,1:end-1)=Nq;
    Uq=Lq';
    Nq = Lq(:, 1:end-1) + Uq(:, 2:end);
end

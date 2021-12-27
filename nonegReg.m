function [S,Jout]=nonegReg(Y,A,S,tol)
% Usage: S=nonegReg(Y,A,S)
%
% Solves
%       minimize    1/2*||Y-A*S||^2
%       subject to  S >= 0
%
% Input:
%        Y: n x p data vector
%        A: n x r regression matrix
% Output:
%        S: [rxp] nonnegative regression 
%        Jout: Cost func
% Magnus Ulfarsson, 2012
    ver=0;
    if nargin<4
        tol=1e-4;
    end
    [n,p]=size(Y);
    [n,r]=size(A);
    index=1:r;
    
    if isempty(S)
        S=rand(r,p);
    end
    AtA=A'*A;
    AtY=A'*Y;
    if n<p
        trYtY=trace(Y*Y');
    else
        trYtY=trace(Y'*Y);
    end
    J0=0.5*trYtY-trace(S*AtY')+0.5*trace(S*S'*AtA);
    Jout(1)=J0;
    for i=1:10000
        for j=1:r,
            ind=index;
            ind(j)=[];
            S(j,:)=(AtY(j,:)-AtA(j,ind)*S(ind,:))/AtA(j,j);   
            S(j,:)=S(j,:).*(S(j,:)>0);            
        end
        Jout(i+1)=0.5*trYtY-trace(S*AtY')+0.5*trace(S*S'*AtA);
        if (Jout(end-1)-Jout(end))/Jout(end) < tol
            break;
        end
        %if rem(i,1)==0
        %    Rj=Y-A*S;
        %    hist(Rj(:),100);shg
        %end
    end
    if(ver==1)
        disp(['Number of iterations: ' num2str(i)]);
        disp(['First order conditions: ']);
        A'*(y-A*s)
    end
end
   

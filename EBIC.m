function EBIC=EBIC(Y,S,A,alpha)
%
%
[T,~]=size(Y);
[M,r]=size(A);
sigma2=1/(T*M)*norm(Y-S*A','fro')^2;
de=sum(S(:)>0)+M*r-r^2;
EBIC=M*log(sigma2)+M+1/T*(log(T)+4*alpha*log(M))*de;
end
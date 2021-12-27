function [B,I,J]=updateB(A,B,k)
% Creates bilinear endmember combos using Mxr endmember matrix 
% A, where M are the number spectral bands and r is the number of endm.
    [~,r]=size(A);
    N=r*(r-1)/2;
    I=zeros(1,N);
    J=I;
    idx=0;
    for i=1:r-1
        for j=i+1:r
            idx=idx+1;
            if j==k || i==k
                B(:,idx)=A(:,i).*A(:,j);
                I(idx)=i;J(idx)=j;
            end
        end
    end
end
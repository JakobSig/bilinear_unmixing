function S=sparseBilinearUnmixingSStep(Y,Ad,Bd,S,hS,hZ,q,stopCriteria,checkConv)
    
    stopCriteria=stopCriteria/10;
    A=[Ad Bd];
    H=[hS hZ];
    H=diag(H)*q;
    ATA=A'*A;
    YA=Y*A;
    for Siter=1:10000
        Sold=S;
        S=(S.*YA)./(S*ATA+(max(eps,S).^(q-1)*H));
        if rem(Siter,checkConv)==0
            ndS=norm(S-Sold,'fro')/norm(Sold,'fro');
            if(ndS<stopCriteria)
                return;
            end
        end
    end
end

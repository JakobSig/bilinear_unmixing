function J=calcBilFanCost(Y,S,A,G,q,h,delta,u)
% J= 1/2||Y-SA'-ZB'||^2_F+h*||S||^q + delta^2/2 ||1'-1'S||+u/2||G||^2_F
    B=createB(A);
    Z=createB(S);
    [P,M]=size(Y);
    [M,r]=size(A);
    [M,N]=size(B);
    X=[Y delta*ones(P,1)]-S*[A; delta*ones(1,r)]'-(G.*Z)*[B;zeros(1,N)]';
    J=0.5*sum(X(:).^2);
    if h>0
        J=J+h*sum(S(:).^q);
    end
    %if delta>0
    %    J=J+0.5*delta^2*sum((sum(S,2)-1).^2);
    %end
    if u>0
        J=J+0.5*u*sum(G(:).^2);
    end
end

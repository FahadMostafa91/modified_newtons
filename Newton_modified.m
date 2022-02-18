function [x,k] = Newton_modified(f,Hessian,x0,tol,kmax) 
x=x0;
n=length(x);
k = 1;
c=0.5;
%% modified  Newton's method
m = f(x); % m=f(xk)

s = norm(m)^2; % s= ||f(xk)||^2
%% estimate hessian matrix 
while (s > tol) && (k < kmax)
    
  Hk = Hessian(x);
    
    [Q,R] = qr(Hk); 
    dk = -backward(R,Q'*m);
    
%    s = m'*m;       
       
    Dh = m'*Hk;        
    mj = f(x+dk);  

    % solve for alpha of 2^j
    L = mj'*mj;                     
    z = c*norm(dk)*norm(Dh); 
    Rh = s - z;  
    j = 0;
    Lmin = L; 
    index = 0;

    while (L > Rh) && (j < 20)
        j=j+1;
        mj = f(x+2^(-j)*dk); 

        L = mj'*mj;
        
        Rh = s - 2^(-j)*z;
        if L < Lmin
            Lmin = L;
            index = j;
        end
    end
    x = x + 2^(-index)*dk;    
    k = k + 1;
    m = f(x);

    s =  norm(m)^2;
end

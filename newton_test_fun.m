clear;
x = [10;10];
xlength = length(x);
a = 2; 
b = 100;
f = @(x) b*(x(2)-x(1)^2)^2+(x(1)-a)^2;
Gradient = @(x) [-4*b*(x(2)-x(1)^2)*x(1)+2*(x(1)-a);
                2*b*(x(2)-x(1)^2)];
Hessian = @(x) [-4*b*(x(2)-x(1)^2)+8*b*x(1)^2+2,-4*b*x(1);
                 -4*b*x(1),2*b];

count = 0;
kmax = 100000; tol = 1e-6; n= length(x); epsilon = 1e-6;

[x,k] = Newton_modified(Gradient,Hessian,x,tol,kmax);
disp(x)

close all
N = 51;
x = linspace(-1.2,1.2,N);
y = x;
z = zeros(N,N);
for i= 1:N
     for j=1:N
         z(i,j) = 100*(y(j)-x(i)^2)^2+(a-x(i))^2;
     end
 end
 surf(x,y,z)
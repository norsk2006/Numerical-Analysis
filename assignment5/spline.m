%%%%%%%%%%%%%%%%%%%%%    compute spline coefficients here    %%%%%%%%%%%%%%%

format long e

%compute a vector of 100 evenly-spaced points between 0 and 2pi, including the endpoints
x = linspace(0, 2*pi,100);

%evaluate vector of  "exact" values of cos at above x values; use MATLAB's built-in cos function
y = arrayfun(@(x) cos(x), x);


%compute the values in vector b here

n=100;


rhs(1) = 0;
for i=2:n-1
    delta = x(i) - x(i-1);
    rhs(i) = (2/delta)*(cos(x(i))-cos(x(i-1)));
end





%compute vector of values in vector a here

A_matrix = zeros(n-1, n-1);
A_matrix(1, 1) = 1;
for i=2:n-1
        A_matrix(i, i) = 1;
        A_matrix(i, i-1) = 1;
   
end
size(A_matrix)

for i=1:n-1
    a(i) = cos(x(i));
end

%compute vector of values in c here 
rhs = rhs';
size(rhs)
size(A_matrix)
b = forsub(A_matrix, rhs)
b = b';
for i=1:n-1
    delta = x(i+1) - x(i);
    c(i) = (cos(x(i+1))-cos(x(i)))/(delta^2) - b(i)/delta;

end



%%%%%%%%%%%%%%%%%%call cosine function here - don't modify%%%%%%%%%%%%%%%%%

X1 = 250;

Y1 = cosine(X1,a,b,c)

X2 = -100;

Y2 = cosine(X2,a,b,c)

X3 = 10*rand(1)-100;

Y3 = cosine(X3,a,b,c)

X4 = 10*rand(1)+100;

Y4 = cosine(X4,a,b,c)
%%%%%%%%%%%%%%%%%%%%%plot spline - don't modify%%%%%%%%%%%%%%%%%%%%%

%number of points on which to plot. n = number of nodes
nplot  = (n-1)*19+1;

xplot = zeros(nplot,1);
yplot = zeros(nplot,1);

%spacing between plot points
nspace = (x(n)-x(1))/(nplot-1);

k = 0;
for i = 1:n-1
    
    for j = 1:19
        
        k = k+1;
        xplot(k) = x(i) + (j-1)*nspace;
        yplot(k) = a(i) + b(i)*(xplot(k) - x(i)) + c(i)*(xplot(k) - x(i))^2;
        
    end
    
end

xplot(nplot) = x(n);
yplot(nplot) = a(i) + b(i)*(x(n) - x(n-1)) + c(i)*(x(n) - x(n-1))^2;

plot(xplot,yplot)

figure

abserr = abs(yplot - cos(xplot));

plot(xplot,abserr)


%%%%%%%%%%%% Compute your cosine function here
function y = cosine( X, a, b, c )
%input
%x - where to evaluate sin function
%a, b, and c are vectors containing quadratic spline coefficients

xspace = linspace(0, 2*pi,100);
x0  = mod(X, 2*pi);

if x0 < 0
    x0 = x0*(-1);
end

for i=1:99
    if x0 > xspace(i) & x0 < xspace(i+1)
        k = i;
        break
    end
end

y = a(k) + b(k) * (x0-xspace(k)) + c(k) * (x0-xspace(k))^2;

end

%%%%%%%%%%%% forward substitution function
function [y] = forsub(L,b)

%Forward-substitution
%accepts an nX1 vector b, an nXn lower triangular matrix L
%generates an nX1 solution vector y

n = size(b,1);

y = zeros(n,1);

y(1) = b(1);

for i = 2:n 
    
    y(i) = b(i);
    
    for j = 1:i-1 
        
        y(i) = y(i) - L(i,j)*y(j); 
        
    end 
    
end

end
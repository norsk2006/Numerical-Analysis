%Example 5.5.1
A = [1 1 3; 
    0 2 1;
    0 0 1;
    -1 -1 -1];

StdGS(A)
ModGS_A = mgson(A)

condition(A)



function U = ModGS(X)
    
    [~, m] = size(X);
    
    U(:,1) = X(:,1)/norm(X(:,1));
    
    U(:,2:m) = X(:,2:m);
    for k = 2:m
            
            
            for j =k:m   
            U(:,j) = U(:,j) - (ctranspose(U(:,k-1))*U(:,j))*U(:,k-1);

            
            
            end
            U(:,k) = U(:,k)/norm(U(:,k));
            
    end

end


function [Q, R] = mgson(X)

[d,n] = size(X);
m = min(d,n);
R = zeros(m,n);
Q = zeros(d,m);
for i = 1:m
    v = X(:,i);
    for j = 1:i-1
        R(j,i) = Q(:,j)'*v;
        v = v-R(j,i)*Q(:,j);
    end
    R(i,i) = norm(v);
    Q(:,i) = v/R(i,i);
end
R(:,m+1:n) = Q'*X(:,m+1:n);
end

function U = StdGS(X)

    [n, m] = size(X);
    U(:,1) = X(:,1)/norm(X(:,1));

    for i = 2:m
        U(:,i)=X(:,i);
            for j = 1:i-1
                U(:,i) = U(:,i) - (dot(U(:,j), U(:,i))) / (norm(U(:,j)))^2 * U(:,j);
                 U(:,i) = U(:,i)/norm(U(:,i));
            end

           
    end
    

end

function y = condition(Z)

for i = 1:5
    

    k(i) = 2*(i-1);
    A= [1 1 1; 
        10^(-k(i)) 10^(-k(i)) 0; 
        10^(-k(i)) 0 10^(-k(i))];

    B  = [1 1 1;
         1+10^(-k(i)) 1 1;
         1 1+10^(-k(i)) 1;
         1 1 1+10^(-k(i))];
    
    condA(i) = cond(A, 'fro');
    condB(i) = cond(B, 2);

    StdGS_A = StdGS(A);
    StdGS_B = StdGS(B);

    normE_StdGS_A(i) = norm(transpose(StdGS_A)*StdGS_A - eye(3), 'fro');
    normE_StdGS_B(i) = norm(transpose(StdGS_B)*StdGS_B - eye(3), 'fro');

    ModGS_A = ModGS(A);
    ModGS_B = ModGS(B);
    
    E = (ModGS_A)'*ModGS_A;
    E = E - eye(3);
    normE_ModGS_A(i) = norm(E, 'fro');
    E = (ModGS_B)'*ModGS_B;
    E = E - eye(3);
    normE_ModGS_B(i) = norm(E, 'fro');

    
    

end
output = [k; condA; condB; normE_StdGS_A; normE_StdGS_B; normE_ModGS_A; normE_ModGS_B];


sprintf('%2d \t %0.1e \t %.1e \t %0.1e \t %0.1e \t %0.1e \t %.1e\n',output)

end
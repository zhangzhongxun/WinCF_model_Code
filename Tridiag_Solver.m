% Function to solve the tridiagonal linear system

% Input:
% d -> the main diagonal of the coefficient matrix
% e -> the first subdiagonal of the coefficient matrix
% f -> the first superdiagonal of the coefficient matrix
% x -> the right hand side of the equations

% Output:
% y -> solution of the linear system

function y = Tridiag_Solver(d,e,f,x)


n = length(x); % n is the number of unknowns

b=zeros(n,1);
g=zeros(n,1);
a=zeros(n-1,1);

b(1)=d(1);
for j=2:n
   a(j-1)=e(j-1)/b(j-1);
   b(j)=d(j)-a(j-1)*f(j-1);
end

%g(1)=x(1);
for j=2:n
      x(j)=x(j)-a(j-1)*x(j-1);
end

x(n)=x(n)/b(n);
for j=n-1:-1:1
   x(j)=(x(j)-f(j)*x(j+1))/b(j);          
end

y=x(:);  

end


   
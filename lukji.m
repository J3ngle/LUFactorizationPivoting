
%Row oriented LU Book Page 79 program 4 modified
function[A]=lukji(A)
%Takes our matrix from Driver
[n,m]=size(A);

%Stops our code if A is not a square matrix
if n~=m, error ('Only square systems');
end
%Perform LU by rows
for k=1:n-1
    if A(k,k)==0
        error('Null pivot element');
    end
    A(k+1:n,k)=A(k+1:n,k)/A(k,k);
    for j=k+1:n
        i= k+1:n;
        A(i,j)=A(i,j)-A(i,k)*A(k,j);
    end
end
return


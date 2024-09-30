%diflu book page 79
%Program 5
%Lu factorization of matrix A jki version
function[A]=diflu(A)
%LUJKI LU factorization of a matrix A in the jki version
[n,m]=size(A);
if n~=m
    error('Only square systems');
end
for j=1:n
    if A(j,j)==0
        error('Null pivot Element');
    end
    for k=1:j-1
        i=k+1:n;
        A(i,j)=A(i,j)-A(i,k)*A(k,j);
    end
    i=j+1:n;
    A(i,j)=A(i,j)/A(j,j);
end
return

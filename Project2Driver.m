%% Main driver code
%Jonathan Engle
%Project 2
%Due 10/23/23
clear;clc;clear all 
tic
%% Step 1
%Generate a random A
global n ;
%n=20
%n=200
%n=1000
n=200;
%trial used 3 
%n=3;
A=randi([-100,100],[n,n]);
checka=rcond(A)>eps;

%Create symetric positive definite matrix
lil=tril(randi([1,1000],[n,n]));
Spd=lil*transpose(lil);
%A=lil;
%To run the symetric positive definite matrix,
%set A to lil, uncomment bellow to check
%A=lil;

if checka==0;
    error('Singular Matrix, try again');
end

if det(A)==0;
    error('Singular Matrix, try again')
end
%This calls our function and tells us which method to proceed with
decision(A)
%pause

%% Step 2
%%Step 2 Given a matrix A, we need to generate a vector b basaed
%on a chosen solution x. I.E for a given A choose a random x and
%generate b via A*x

%Generate random vector x
rng("default")
x = rand(n,1);
%compute b
global b;
b=A*x;

%% Step 3
% Test Matrix for :
%No pivot, partial pivot, Complete pivot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Given Matrix tested
Atest=[3 17 10; 2 4 -2; 6 18 -12];
%No pivoting using modified book algorithm
lukji(Atest);
%Partial pivoting
partialPivoting(Atest);
[L, U, P] = partialPivoting(A);
%Generating a blank matrix
storp=eye(n);
storp=storp+U-2*eye(n);
%Storing our L and U solutions collectively to make it More efficient
storp=L+storp;
%Complete pivot
completePivoting2(Atest);
% storc=eye(n);
% storc=storc+Ucomp-2*eye(n);
% storc=Lcomp+storc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%For a randomly generated matrix
%No pivot
nopiv=spiv(A);
[spL, spU] = spiv(A);

%Partial Pivot
[L, U, P] = partialPivoting(A);
ppA=partialPivoting(A);

%Complete pivot
[Lcomp, Ucomp, Pcomp, Qcomp] = completePivoting2(A);
completePivoting2(A);
storreal=eye(n);
storreal=storreal- 2*eye(n);

%% Step 4
%Forward substitution for our random matrix
y=forwardSubstitution(spL,b);
%Backward substitution for our random matrix
xnew=backwardSubstitution(spU,y);

%%% Uncomment when you need to use partial Pivoting
%Forward substitution for our random matrix
ypp=forwardSubstitution(L,b);
%Backward substitution for our random matrix
xpp=backwardSubstitution(U,ypp);

%%% Uncomment when you need to use partial Pivoting
%Forward substitution for our random matrix
ycc=forwardSubstitution(Lcomp,b);
%Backward substitution for our random matrix
xcc=backwardSubstitution(Ucomp,ycc);


%% Step 5
%Checking the accuracy

%Errors for no pivoting
norm(A-spL*spU,1)/norm(A,1);
norm(A-spL*spU,"fro")/norm(A,"fro");
%Errors for partial pivoting
norm(P*A-L*U,1)/norm(A,1);
norm(P*A-L*U,"fro")/norm(A,"fro");
%Errors for Complete pivoting
norm(Pcomp*A*Qcomp-Lcomp*Ucomp,1)/norm(A,1);
norm(Pcomp*A*Qcomp-Lcomp*Ucomp,"fro")/norm(A,"fro");

%Accuracy check for x
norm(x-xnew,1)./norm(x,1); %One norm error
norm(x-xnew,2)./norm(x,2); %2 Norm error

%Residual error for no pivoting
norm(b-A*xnew,1)/norm(b,1);
norm(b-A*xnew,2)/norm(b,2);

%Residual error for partial pivoting
norm(b-P*A*xpp,1)/norm(b,1);
norm(b-P*A*xpp,2)/norm(b,2);

%Residual error for complete pivoting
norm(b-Pcomp*A*Qcomp*xcc,1)/norm(b,1);
norm(b-A*xcc,2)/norm(b,2);
%% Correctness Test Task
%Given the matrix and vecot
Acorrect= [2 1 0;-4 0 4; 2 5 -10];
bcorrect= [3 ;0 ;17];

%Partial Pivoting 
[Lco, Uco, Pco] = partialPivoting(Acorrect);
ppAcorrect=partialPivoting(Acorrect);
%Forward substitution for our given matrix
ycomp=forwardSubstitution(Lco,bcorrect);
%Backward substitution for our given matrix
xcomp=backwardSubstitution(Uco,ycomp)
%Above solves for our x
%This part is just a simple check that our algorithm works
o=Acorrect*xcomp;
Pco*o;
Lco*Uco*xcomp-bcorrect;
%Woohoo!
decision(Acorrect);
%%Complete pivoting 
[Lcomp2, Ucomp2, Pcomp2, Qcomp2] = completePivoting2(Acorrect);
completePivoting2(Acorrect);
%Forward substitution for our given matrix
ycomp2=forwardSubstitution(Lcomp2,bcorrect);
%Backward substitution for our given matrix
xcomp2=backwardSubstitution(Ucomp2,ycomp2)
Lcomp2*Ucomp2*xcomp
Qcomp2;

toc
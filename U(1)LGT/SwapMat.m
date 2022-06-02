function P=SwapMat(m,n)
    %----------------------------------------------------------------------
    % Let A and B be two Hilbert spaces of dimensions m and n,
    % respectively. The SwapMat is a linear map from A\otimes B to B\otimes
    % A such that |i>|j> is mapped to |j>|i>. The basis for A\otimes B or
    % B\otimes A are chosen according to the canonical Kronecker product
    % convention. For example, basis vectors of A\otimes B are ordered as 
    % |1>|1>,|1>|2>,...,|1>|n>,|2>|1>,|2>|2>,...,|2>|n>,...,|m>|n>. 
    %----------------------------------------------------------------------

    Dim=m*n; % Total Hilbert space dimension. 
    P=sparse(1:Dim,reshape(reshape(1:Dim,[n,m]).',[1,Dim]),1,Dim,Dim);
end
%{
% Test. 
P=SwapMat(2,3);
u1=kron([1.1;2.2],[pi;-3;5*1i]);
v1=kron([pi;-3;5*1i],[1.1;2.2]);
disp(v1-P*u1);
uA2=rand(2,1);
uB2=rand(3,1);
u2=kron(uA2,uB2);
v2=kron(uB2,uA2);
disp(v2-P*u2);
%}
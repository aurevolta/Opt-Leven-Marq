function [r,a_new] = OLM_ROBUST_STEP_2(parameters,B,U1,g,D,lambda)
%OLM_ROBUST_STEP_2 Second step of the robust LM step
%
% [r,a_new] = LM_ROBUST_STEP_2(a,na,data,n,B,U1,g,D,lambda) computes the
% residual [r] and the new state [a_new] using the output
% from the step 1 [B], [U1], [g], [D], [lambda] as well as the state [a],
% its dimension [na], the data set [data] and the number of
% equations/measurements. 
%
% see also OLM_ROBUST_STEP_1, OLM_ROBUST_STEP_LAMBDA

% SPDX-License-Identifier: Apache-2.0
% 2016 Aureliano Rivolta

%%

% QR decomposition of the transformed problem
[Q2,U2] = qr([U1;B*lambda]);

% compute u
Q2T = Q2';
u = Q2T(1:parameters.na,1:parameters.n)*g;

% get the solution (backward substitution) to b2=-U2(1:na,1:na)^-1*u;
b2 = backward_subs(U2(1:parameters.na,1:parameters.na),-u,parameters.na);

% descale the output
delta_a = D*b2;

a_new =parameters.a+delta_a;

% compute the residual of the found solution (but not the jacobian)
[r] =parameters.fun(a_new,parameters.data,1);

end

function  x = backward_subs(U,b,n)
% BACKWARD_SUBS solve upper triangular problem
%
% x = BACKWARD_SUBS(U,b,n) computes [x], the solution of U * x = b with [U]
% upper triangular matrix, the vector [b] and the size of the problem [n]. 

% Author: Aureliano Rivolta
% e-mail: aureliano.rivolta@polimi.it
% year: 2016
% private use only

%%

% n is the size of the problem
x = zeros(n,1);

% compute the last entry of the solution vector x
x(n) = b(n)/U(n,n);

% substitute from the last-1 entry backward to the first
for i = (n-1) : -1 : 1
    x(i) = (b(i)-sum(U(i,i:n)*x(i:n)))/U(i,i);
end

end
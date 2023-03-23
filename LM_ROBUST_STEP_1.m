function [B,U1,g,D,Q1T,r] = LM_ROBUST_STEP_1(parameters)
% LM_ROBUST_STEP_1 First step of robust LM
%
% [B,U1,g,D,Q1T] = LM_ROBUST_STEP_1(a,na,data,n) computes the matrices [B],
% [U1] and [Q1T] and vector [g] of transformed residuals for the LM step 2
% including also the Jacobian scaling [D] from the stored parameters 
% structure [parameters] 
%
% see also LM_ROBUST_STEP_2, LM_ROBUST_STEP_LAMBDA

% SPDX-License-Identifier: Apache-2.0
% 2016 Aureliano Rivolta

%%
[r,J] = parameters.fun(parameters.a,parameters.data,1);

% normalize the Jacobian
[Jscaled,D] = scaleJacobian(J);

% compute the QR decomposition of the Jacobian
[Q1,U1] = qr(Jscaled);

% compute the forcing term from residual
Q1T = Q1';
g = Q1T*r;

% compute the second term (as result of scaling)
% A=Jscaled'*Jscaled;%no scaling
% B=diag(sqrt(diag(A))); %no scaling
B = eye(parameters.na);

end
function [a_new,r,lambda] = OLM_FAST_STEP(parameters)
% OLM_FAST_STEP Fastest version of LM single step
%
% [a_new,r,lambda] = OLM_FAST_STEP(a,na,data,n,compute_r,mu0,W,r_old)
% computes the new iteration [a_new] with the new residual [r] and the
% updated damping [lambda] from the current state estimate [a], its
% dimensions [na], the data seta [data], the number of
% equations/measurements [n] the flag [compute_r], the initial damping
% [mu0] and the previous residual [r_old]
% 
% see also OLM_FAST, OLM_SCALEJACOBIAN

% SPDX-License-Identifier: Apache-2.0
% 2016 Aureliano Rivolta

%%


% compute residual and Jacobian if needed
if parameters.compute_r
    [ru,J] = parameters.fun(parameters.a,parameters.data,parameters.compute_r);
    r = ru(1:parameters.n,1);
else
    [~,J] = parameters.fun(parameters.a,parameters.data,0);
    r = parameters.r(1:parameters.n,1); % take last computed residual as correct for the computation
end

% determine the damping parameter
lambda = (parameters.mu0*norm(r)^parameters.lambda_exponent)/parameters.n;

% scale the Jacobian
[Jscaled,D] = OLM_scaleJacobian(J);

% compute the A matrix
A = Jscaled'*diag(parameters.W)*Jscaled;
% A=A+lambda*diag(diag(A)); % if it is scaled diag(A)=eye(na)
A = A+lambda*eye(parameters.na);

% compute the gradient
b = -Jscaled'*diag(parameters.W)*r;

% computethe solution
delta_a = D*(A\b);

% assemble the step state (1:parameters.na,1) 
a_new= parameters.a(1:parameters.na,1)+delta_a(1:parameters.na,1);

% check the new solution (compute just the residual)
[ru] = parameters.fun(a_new,parameters.data,1);
 r = ru(1:parameters.n,1);
 
end
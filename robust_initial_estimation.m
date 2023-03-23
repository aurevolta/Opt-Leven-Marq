function parameters = robust_initial_estimation(parameters)
% ROBUST_INITIAL_ESTIMATION Estimates the best initial damping parameter
%
% parameters = ROBUST_INITIAL_ESTIMATION(parameters) computes the initial damping parameter
% [mu0] and the weights [W] from the parameter structure [parameters] and
% adjourne the structure.
% [mu0] is computed looking at different initial damping parameters for the
% LM and choosing the one that minimizes the most the residual criterion.
%
% see also LM_ROBUST_STEP_LAMBDA, LM_ROBUST_STEP_1, LM_ROBUST_STEP_2,
% METRIC

% SPDX-License-Identifier: Apache-2.0
% 2016 Aureliano Rivolta

%%

% range of lambda to be tested
LAM = 10.^linspace(-3,1,4);

% test the various lambdas using the robust approach (less computation for
% the same initial conditions and just varying lambda)
[lambda_best,rho_best,r_best] = LM_ROBUST_STEP_LAMBDA(parameters,LAM);

% determine the initial reduction parameter
parameters.mu0 = lambda_best/norm(parameters.r)*parameters.n;

if parameters.removal && rho_best > 0
    % outliers removal if needed
    parameters.W (:,1) = parameters.removal_function(r_best,parameters);
end

end
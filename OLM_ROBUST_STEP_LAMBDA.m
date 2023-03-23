function [lambda_best,rho_best,R] = OLM_ROBUST_STEP_LAMBDA(parameters,lambda)
% OLM_ROBUST_STEP_LAMBDA Routine to compute the best damping step
%
% [lambda_best,rho_best,R] = OLM_ROBUST_STEP_LAMBDA(parameters,lambda)
% computes the best damping [lambda_best] along with the best increment 
% [rho_best] and the relative residual vector [R] from the parameters 
% structure [parameters] 
%
% see also OLM_ROBUST_INITIAL_ESTIMATION, OLM_ROBUST_STEP_1, 
% OLM_ROBUST_STEP_2, OLM_METRIC

% SPDX-License-Identifier: Apache-2.0
% 2016 Aureliano Rivolta


%%

% first step is computed just once
[B,U1,g,D,~] = OLM_ROBUST_STEP_1(parameters);

% initialize 
RHO = zeros(1,length(lambda));
W0 = ones(parameters.n,1);
I = length(lambda);
r = zeros(parameters.n,I);

% computes all the steps
for i = 1 : I
    % second step
    [r(:,i),~] = OLM_ROBUST_STEP_2(parameters,B,U1,g,D,lambda(i));
    % compute the metric
    [~,rho] = OLM_metric(r(:,i),parameters);
    % save the metric
    RHO(i) = rho;
end

% obtain the best metric
rho_best = max(RHO);

% in case of exceptions
if isnan(rho_best)
    lambda_best = 1;
    R = r(:,1);
else
    % find the maximum metric
    INDEX = (rho_best==RHO);
    % find the damping that corresponds to the best metric
    lambda_best = min(lambda(INDEX));
    R = r(:,INDEX);
end

end
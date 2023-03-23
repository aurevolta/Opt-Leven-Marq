function [chi,rho] = OLM_metric(r,parameters)
% OLM_METRIC Metric and advancement computation 
%
% [chi,rho] = OLM_METRIC(r,W,chi_best) computes the metric [chi] and the
% advancement parameter [rho] using the residual vector [r], the weighting
% factors [W] and the best metric [chi_best] achieved in the previous steps
% and stored inside the input structure [parameters]
% Dimensions of [r] and [W] shall be consistent, [chi,rho] are 
% scalars. 
%
% see also OLM, OLM_FAST, OLM_ROBUST

% SPDX-License-Identifier: Apache-2.0
% 2016 Aureliano Rivolta

%%

% chi squared criterion (metric)

% chi=sum(r.^2);
% chi=sum(r.^2)/length(r);
chi = sum((r.*parameters.W).^2)/parameters.n/(sum(parameters.W).^2);
% chi = sum((r.*parameters.W).^2)/(sum(parameters.W));
% chi=sum((r.*W).^2)/length(r);
% chi=sum((r.*W).^2)/length(r)/(sum(W));


% metric function update
rho = (parameters.chi_best-chi);
% rho=0.5*(chi_best-chi)/(delta_a'*(lambda*delta_a+b));


end
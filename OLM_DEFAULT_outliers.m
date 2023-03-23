function  W = OLM_DEFAULT_outliers(r,parameters)
% OLM_DEFAULT_OUTLIERS Default no-removal function
%
% W = DEFAULToutliers_removal_function(r,parameters) computes the
% default weights [W] of the residual [r] from the input parameter strcuture 
% [parameters]
%
% see also 

% SPDX-License-Identifier: Apache-2.0
% 2016 Aureliano Rivolta

%%
W = ones(size(parameters.W));

end
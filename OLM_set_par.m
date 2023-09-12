function parameters = OLM_set_par(a0,fun,version,data,removal_function)
% OLM_SET_PAR Set LM parameter structure
%
% parameters = OLM_SET_PAR(a0,fun,version) construct the structure [parameters] for
% the minimization of the function [fun] with initial guess [a0]. No data,
% equation remove options are considered in this case but the variable
% [version] allows to choose between a robust but slower version ['robust']
% and a computationally faster implementation ['fast'].
%
% parameters = OLM_SET_PAR(a0,fun,version,data) is the same as before but the input
% [data] is added in the structure. No restrictions on the natrue of [data]
%
% parameters = OLM_SET_PAR(a0,fun,version,data,removal_function) adds an outlier
% removal function handle [removal_function] with a precise input-output
% structure.
%
% see also OLM, OLM_FAST, OLM_ROBUST, OLM_DEFAULT_METRIC

% SPDX-License-Identifier: Apache-2.0
% 2016 Aureliano Rivolta

%% input parsing

% states related info
parameters.a = a0;
parameters.na = length(a0);

% add user defined data (if present)
if nargin > 3
    parameters.data = data;
else
    parameters.data = [];
    data = [];
end

% add user defined choice for removal
if nargin > 4
        parameters.removal = 1;
        parameters.removal_function = removal_function;
else
    parameters.removal = 0;
    parameters.removal_function = @OLM_DEFAULT_outliers;
end

% initial computation of user defined function

if isa(fun,'function_handle')
    
    % compute the first residual
    r = fun(a0,data,1);
    
else
    error 'User defined function is not a function'
end

%% initialization

% initialize flags
parameters.compute_r = 1;

% % initialize weights as scalar
% parameters.W = 1;

% initialize LM damping
parameters.lambda = 1;

% store the function
parameters.fun = fun;

% set the LM mode
parameters.version = version;

% set the default criterion
parameters.metric = @OLM_DEFAULT_metric;

% extract info
parameters.r = r;
parameters.n = length(r);
parameters.W = ones(length(r),1);

%% Default parameters

% set the maximum number of iterations
parameters.n_iter = 50;

% decreasing LM damping law parameter exponent
parameters.lambda_exponent = 1.05;

% initial lambda scaling
% parameters.mu0 = 1;
parameters.mu0 = 1/(r'*r);

% maximum consecutive failed steps
parameters.max_failures = 3;

% stoppinc criteria
parameters.chi_threshold = 1e-24;

% damping for outliers removal
parameters.damping = 1;

%% computes the initial chi squared criterion
parameters.chi_best = zeros(1,'like',r);

% parameters.chi_best = OLM_DEFAULT_metric(r,parameters);
parameters.chi_best = parameters.metric(r,parameters);

end







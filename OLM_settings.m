function parameters = OLM_settings(varargin)
% OLM_SETTINGS Set LM parameter structure (NO CODEGEN)
%
% parameters = OLM_SETTINGS(a0,fun,version) construct the structure [parameters] for
% the minimization of the function [fun] with initial guess [a0]. No data,
% equation remove options are considered in this case but the variable
% [version] allows to choose between a robust but slower version ['robust']
% and a computationally faster implementation ['fast'].
%
% parameters = OLM_SETTINGS(a0,fun,version,data) is the same as before but the input
% [data] is added in the structure. No restrictions on the natrue of [data]
%
% parameters = OLM_SETTINGS(a0,fun,version,data,removal_function) adds an outlier
% removal function handle [removal_function] with a precise input-output
% structure.
%
% see also OLM, OLM_FAST, OLM_ROBUST, OLM_DEFAULT_METRIC, OLM_SET_PAR

% SPDX-License-Identifier: Apache-2.0
% 2016 Aureliano Rivolta

%% input parsing

% generates parser object
p = inputParser;

% add required parameters
addRequired(p,'a0');
addRequired(p,'fun',@(x) isa(x,'function_handle'));
addRequired(p,'version',@(x) ischar(x) && ismember(x,{'fast','robust'}));

% add optional parameters or functions
addParameter(p,'data',[]);
addParameter(p,'removal_function',@OLM_DEFAULT_outliers,...
    @(x) isa(x,'function_handle') && nargin(x)==2);
addParameter(p,'metric',@OLM_DEFAULT_metric,...
    @(x) isa(x,'function_handle') && nargin(x)==2);

% parse the user input
parse(p,varargin{:})

%% extract user defined parameters

% first guess of the solver
parameters.a = p.Results.a0;

% additional data for the function evaluation
parameters.data = p.Results.data;

% function to be minimized
parameters.fun = p.Results.fun;

% removal function
parameters.removal_function = p.Results.removal_function;

% set the LM mode
parameters.version = p.Results.version;

% set the LM chi squared criterion
parameters.metric = p.Results.metric;

% states related info
parameters.na = length(parameters.a);

% compute the first residual (with jacobian)
[r,J] = parameters.fun(parameters.a,parameters.data,1);

% residual save and number of equations in the minimizer
parameters.r = r;
parameters.n = length(r);

% check sizes of the jacobian
[J_n,J_na] = size(J);
if J_n ~= parameters.n || J_na ~= parameters.na
    error('Jacobian size is [%.0f %.0f], expected [%.0f %.0f]',...
        J_n,J_na,parameters.n,parameters.na);
end

% Weights initialization
parameters.W = ones(parameters.n,1);

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

% initialize residual computation flags
parameters.compute_r = 1;

% initialize LM damping
parameters.lambda = 1;

% by default does not apply removal
parameters.removal = 0;

%% computes the initial chi squared criterion
parameters.chi_best = zeros(1,'like',r);

% parameters.chi_best = OLM_DEFAULT_metric(r,parameters);
parameters.chi_best = parameters.metric(r,parameters);

end







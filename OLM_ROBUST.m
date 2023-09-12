function [a_best,chi_best,W,n_iterations,CHI,A,RESIDUAL,LAMBDA] = OLM_ROBUST(parameters)
% OLM_ROBUST LM algorithm robust version
%
% [a_best,chi_best,W,n_iterations,CHI,A,RESIDUAL,LAMBDA] = OLM_FAST(parameters) 
% computes the best estimate of the state [a_best], its chi squared criterion 
% [chi_best] as well as the weights [W], the number of iterations 
% [n_iterations], the story of criterion [CHI], state [A], residual
% [RESIDUAL] and damping paramters [LAMBDA]. The latest are matrices of 
% dimensions [n] times [n_iter] but carries info until [n_iter]. 
% It requires as input the structure of paramters [parameters] created with
% the corresponding function OLM_SET_PAR
%
% see also OLM_FAST_STEP, OLM_DEFAULT_METRIC, OLM_SET_PAR

% SPDX-License-Identifier: Apache-2.0
% 2016 Aureliano Rivolta

%%

% force residual computation if necessary
parameters.compute_r = 1; 

% hystory values recording variables
CHI = zeros(1,parameters.n_iter);
LAMBDA = zeros(1,parameters.n_iter);
A = zeros(parameters.na,parameters.n_iter);
RESIDUAL = zeros(parameters.n,parameters.n_iter);

% initialize the counters
counter = 0;
compute1 = 1;
iter_counter=0;

%initialize
[B,U1,g,D] = OLM_ROBUST_STEP_1(parameters);

% actual computation
for i=1:parameters.n_iter
    
    % adjourne the iteration counter
    iter_counter = i;    
    
    % two step
    if compute1
        [B,U1,g,D,~,parameters.r] = OLM_ROBUST_STEP_1(parameters);
    end
    
    lambda = (parameters.mu0*norm(parameters.r)^parameters.lambda_exponent)/parameters.n;
    
%     lambda = (parameters.mu0); % ADD
    
    [r_new,a_new] = OLM_ROBUST_STEP_2(parameters,B,U1,g,D,lambda);

    
    % compute the metric
    [chi,rho] = parameters.metric(r_new,parameters);
    
    % save data for diagnostics
    CHI(i) = chi; 
    A(1:parameters.na,i) = a_new;
    RESIDUAL(1:parameters.n,i) = r_new;
    LAMBDA(i) = lambda;
    
    % check if the step is accepted or not
    if rho > 0
        % step accepted
        parameters.chi_best = chi;
        parameters.a(1:parameters.na,1) = a_new(1:parameters.na,1);
        parameters.r(1:parameters.n,1) = r_new(1:parameters.n,1);
        
        % outliers removal
        if parameters.removal
            parameters.W = parameters.removal_function(r_new,parameters);
        end
        
        % reset the failure counter to zero
        counter = 0;
        
        % set compute_r to 0 to reduce recomputation
        parameters.compute_r = 0;
        compute1 = 1;
        
        parameters.mu0 = parameters.mu0/2;
        
        % exit when under threshold
        if parameters.chi_best <= parameters.chi_threshold 
            break
        end
        
        
    else
        % step rejected
        counter = counter+1;
        compute1 = 0;
        parameters.mu0 = parameters.mu0*2;

        if counter == parameters.max_failures
            % terminates after n consecutive missed steps
            break
        end

    end
    
end

%% outputs 
a_best = parameters.a;
chi_best = parameters.chi_best;
W = parameters.W;

% save the number of iterations
n_iterations = iter_counter;

end

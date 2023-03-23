function [a_best,chi_best,W,n_iterations,CHI,A,RESIDUAL,LAMBDA] = LM_FAST(parameters)
% LM_FAST LM algorithm fastest version
%
% [a_best,chi_best,W,n_iterations,CHI,A,RESIDUAL,LAMBDA] = LM_FAST(parameters) 
% computes the best estimate of the state [a_best], its chi squared criterion 
% [chi_best] as well as the weights [W], the number of iterations 
% [n_iterations], the story of criterion [CHI], state [A], residual
% [RESIDUAL] and damping paramters [LAMBDA]. The latest are matrices of 
% dimensions [n] times [n_iter] but carries info until [n_iter]. 
% It requires as input the structure of paramters [parameters] created with
% the corresponding function SET_LM_PAR
%
% see also LM_FAST_STEP, METRIC, SET_LM_PAR

% SPDX-License-Identifier: Apache-2.0
% 2016 Aureliano Rivolta

%%


parameters.compute_r = 1; % force if necessary

% hystory values recording variables
CHI = zeros(1,parameters.n_iter);
LAMBDA = zeros(1,parameters.n_iter);
A = zeros(parameters.na,parameters.n_iter);
RESIDUAL = zeros(parameters.n,parameters.n_iter);

% initialize the counters
counter=0;
iter_counter=0;

% actual computation
for i=1:parameters.n_iter
    
    % adjourne the iteration counter
    iter_counter = i;
    
    % use the fast LM step
    [a_new,r_new,lambda] = LM_FAST_STEP(parameters);
    
    % compute the metric
    [chi,rho] = metric(r_new,parameters);
    
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
        
        parameters.mu0 = parameters.mu0/5; % ADD
        
        % exit when under threshold
        if parameters.chi_best <= parameters.chi_threshold 
            break
        end
        
    else
        % step rejected
        counter = counter+1;
        
        parameters.mu0 = parameters.mu0*5; % ADD

        if counter == parameters.max_failures
            % terminates after 3 consecutive missed steps
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

function [a_best,chi_best,W,n_iterations] = OLM_FAST(parameters)
% OLM_FAST LM algorithm fastest version
%
% [a_best,chi_best,W,n_iterations,CHI,A,RESIDUAL,LAMBDA] = OLM_FAST(parameters) 
% computes the best estimate of the state [a_best], its chi squared criterion 
% [chi_best] as well as the weights [W], and the number of iterations 
% [n_iterations].
% It requires as input the structure of paramters [parameters] created with
% the corresponding function OLM_SET_PAR
%
% see also OLM_FAST_STEP, OLM_DEFAULT_METRIC, OLM_SET_PAR

% SPDX-License-Identifier: Apache-2.0
% 2016 Aureliano Rivolta

%%

% force residual computation if necessary
parameters.compute_r = 1; 

% initialize the counters
counter      = 0;
iter_counter = 0;

% actual computation
for i=1:parameters.n_iter
    
    % adjourne the iteration counter
    iter_counter = i;
    
    % use the fast LM step
    [a_new,r_new,~] = OLM_FAST_STEP(parameters);
    
    % compute the metric
    [chi,rho] = parameters.metric(r_new,parameters);
    
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
        counter = counter + 1;
        
        parameters.mu0 = parameters.mu0*5; % ADD

        if counter == parameters.max_failures
            % terminates after 3 consecutive missed steps
            break
        end

    end
    
end

%% outputs 

% save solution, error and weights
a_best          = parameters.a;
chi_best        = parameters.chi_best;
W               = parameters.W;

% save the number of iterations
n_iterations    = iter_counter;

end

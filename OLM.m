function [varargout] = OLM(parameters,FGS)
% OLM Levenberg-Marquardt minimization algorithm
%
% [a_best,chi_best,W,n_iterations,CHI,A,RESIDUAL,LAMBDA] = OLM(parameters)
% computes the best estimate of the state [a_best], its chi squared criterion
% [chi_best] as well as the weights [W], and the number of iterations
% [n_iterations]. If the robust version is selected it will also output the
% story of criterion [CHI], state [A], residual [RESIDUAL] and damping
% paramters [LAMBDA]. The [parameters] structure contains all the info on
% the problem and the Levenberg-Marquardt algorithm options, including also
% the version: 'fast' or 'robust'. 
%
% [a_best,chi_best,id] = OLM(parameters,FGS) calls [n] times the function with
% [n] different initial guesses provided in the new input [FGS]. This
% version picks the best estimate from all the first guesses [FGS] and
% outputs the solution [a_best] and the corresponding fitting criterion
% [chi_best]. Others parameters are not computed. A robust initial
% estimation of parameters is performed. [id] returns the id of the best
% estimate.
%
%
%
% see also OLM_FAST, OLM_ROBUST, OLM_SET_PAR

% SPDX-License-Identifier: Apache-2.0
% 2016 Aureliano Rivolta


%%

if nargin == 1
    
    % switch between modes
    if strcmp(parameters.version,'fast')
        
        % call the fast version
        [a_best,chi_best,W,n_iterations] = OLM_FAST(parameters);
        
        % output handling
        varargout{1} = a_best;
        varargout{2} = chi_best;
        varargout{3} = W;
        varargout{4} = n_iterations;
        
    elseif strcmp(parameters.version,'robust')
        
        % call the robust version
        [a_best,chi_best,W,n_iterations,CHI,A,RESIDUAL,LAMBDA] = OLM_ROBUST(parameters);
        
        % eliminates the extra entries in the iteration logs
        CHI         = CHI(1:n_iterations);
        A           = A(:,1:n_iterations);
        RESIDUAL    = RESIDUAL(:,1:n_iterations);
        LAMBDA      = LAMBDA(1:n_iterations);
        
        % output handling
        varargout{1} = a_best;
        varargout{2} = chi_best;
        varargout{3} = W;
        varargout{4} = n_iterations;
        varargout{5} = CHI;
        varargout{6} = A;
        varargout{7} = RESIDUAL;
        varargout{8} = LAMBDA;
        
    else
        error 'wrong version selected'
    end
    
else
    % multiple trials, pick the best one.
    
    % number of trials
    n = size(FGS,2);
    
    % initialize
    a = zeros(size(parameters.a,1),n);
    chi = zeros(1,n);
    
    for i = 1 : n
        % load first guess
        parameters.a = FGS(:,i);
        
        % set robust initial guess
        parameters = OLM_robust_initial_estimation(parameters);
        
        % compute again
        if strcmp(parameters.version,'fast')
            [a(:,i),chi(i)] = OLM_FAST(parameters);
        elseif strcmp(parameters.version,'robust')
            [a(:,i),chi(i)] = OLM_ROBUST(parameters);
        end
        
    end
    
    % find the best estimate
    [~,index]   = min(log10(chi));
    chi_best    = chi(index(1));
    a_best      = a(:,index(1));
    
    % output the id
    id = index(1);
    
    % handle function output
    varargout{1} = a_best;
    varargout{2} = chi_best;
    varargout{3} = id;
    
end


end


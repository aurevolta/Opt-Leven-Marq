clc,close all,clear all


%% test 01

a0 = rand*100;
fun = @testfun_A;
version = 'fast';

tic
parameters = OLM_set_par(a0,fun,version);


parameters = OLM_robust_initial_estimation(parameters);


[a_best,chi_best,W,n_iterations,CHI,A,RESIDUAL,LAMBDA] = OLM(parameters);
toc
figure,semilogy(CHI)
a_best-1
n_iterations


tic
parameters = OLM_set_par(a0,fun,version);
parameters = OLM_robust_initial_estimation(parameters);

parameters.n_iter=300;
parameters.max_failures=10;
parameters.version = 'robust';

[a_best,chi_best,W,n_iterations,CHI,A,RESIDUAL,LAMBDA] = OLM(parameters);
toc
figure,semilogy(CHI)
a_best-1
n_iterations

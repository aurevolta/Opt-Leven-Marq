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



%%

fun = @testfun_B;
version = 'fast';
parameters = OLM_set_par(rand,fun,version);

FGS = linspace(-4*pi,4*pi,9)+1e-12;
% FGS = linspace(pi/2,4*pi,5);


[a_best,chi_best,~,~,~,~,~,~,id] = OLM(parameters,FGS);

fcn = @(a) -sin(a)./a + 0.5;

figure
plot(linspace(-4*pi,4*pi,50),fcn(linspace(-4*pi,4*pi,50)))
hold on
plot(FGS,fcn(FGS),'ko')
plot(FGS(id),fcn(FGS(id)),'*')

plot(a_best,testfun_B(a_best,[],1),'rx')



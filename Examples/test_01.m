clc,close all,clear all


%% test 01

a0 = rand*100;
fun = @testfun_A;
version = 'fast';

tic
parameters = OLM_set_par(a0,fun,version);
parameters = OLM_robust_initial_estimation(parameters);
[a_best,chi_best,W,n_iterations] = OLM(parameters);
toc

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


[a_best,chi_best,id] = OLM(parameters,FGS);

fcn = @(a) -sin(a)./a + 0.5;

figure
plot(linspace(-4*pi,4*pi,100),fcn(linspace(-4*pi,4*pi,100)))
hold on
plot(FGS,fcn(FGS),'ko')
plot(FGS(id),fcn(FGS(id)),'*')
plot(a_best,testfun_B(a_best,[],1),'rx')
legend('function','Initial guesses','best starting point','solution','location','best')
grid on

%% time comparison

fun = @(x) testfun_A(x,[],1);
version = 'fast';
opt = optimoptions('fsolve','algorithm','levenberg-marquardt',...
    'display','none','Jacobian','on','tolx',1e-14,'FunctionTolerance',1e-14);

[T1,T2,T3,E1,E2,E3]=deal(zeros(1,1000));


for i = 1 : length(T1)
    a0 = rand*100;
    
    tic
    sol = fsolve(fun,a0,opt);
    T1(i)=toc;
    E1(i)=abs(1-sol);
    
    tic
    parameters = OLM_set_par(a0,@testfun_A,version);
    parameters = OLM_robust_initial_estimation(parameters);
    [a_best,chi_best,W,n_iterations] = OLM(parameters);
    T2(i)=toc;
    E2(i)=abs(1-a_best);
    
    tic
    parameters = OLM_set_par(a0,@testfun_A,'robust');
    parameters.n_iter=300;
    parameters.max_failures=10;
    parameters = OLM_robust_initial_estimation(parameters);
    [a_best,chi_best,W,n_iterations,CHI,A,RESIDUAL,LAMBDA] = OLM(parameters);
    T3(i)=toc;
    E3(i)=abs(1-a_best);
end

figure
plot(T1,'linewidth',2)
hold on
plot(T2,'linewidth',2)
plot(T3,'linewidth',2)
grid on
axis tight
set(gca,'Xlim',[10,length(T1)])
legend('Matlab fsolve','OLM fast','OLM robust')
% set(gca,'ylim',[0 max([T1,T2,T3])])
title 'Execution time test'


figure
semilogy(E1,'linewidth',2)
hold on
semilogy(E2,'linewidth',2)
semilogy(E3,'linewidth',2)
grid on
axis tight
legend('Matlab fsolve','OLM fast','OLM robust')
% set(gca,'ylim',[0 max([T1,T2,T3])])
title 'Test precision'
ylabel 'error'

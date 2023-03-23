function [r,J] = testfun_A(a,data,compute_r)

g = 1;

if compute_r
    r = (a-g).^2;
else
    r = [];
end

J = 2*(a-g);


end
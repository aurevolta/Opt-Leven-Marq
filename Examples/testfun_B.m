function [r,J] = testfun_B(a,data,compute_r)



if compute_r
    r = -sin(a)./a + 0.5;
else
    r = [];
end

J = (sin(a)-a*cos(a))/(a^2);


end
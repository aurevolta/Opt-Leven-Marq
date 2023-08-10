function [Jscaled,D] = OLM_scaleJacobian(J)
% OLM_SCALEJACOBIAN Jacobian normalization by column
%
% [Jscaled,D] = OLM_SCALEJACOBIAN(J) computes the [n x m]  scaled Jacobian 
% [Jscaled] from the [n x m] Jacobian matrix [J] and the scaling matrix [D]
% such that [Jscaled] = [J]*[D]. Each column of [Jscaled] has unitary norm,
% except for cases where a column is constituted by zeros where scaling is
% not performed. 

% SPDX-License-Identifier: Apache-2.0
% 2016 Aureliano Rivolta

%%

% compute the norm of the i-th column of jacobian
L = sqrt(sum(J.^2,1));

if sum(L==0)>0
    % if one or more column are zero norm (singularity issue) avoid to
    % normalize
    D = eye(length(L));
    Jscaled = J;
else
    % generate the scaling matrix D
    D = diag(L.^-1);
    % scale the Jacobian
    Jscaled = J*D;
end

end
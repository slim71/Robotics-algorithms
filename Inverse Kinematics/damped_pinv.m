function dpinv = damped_pinv(mat, lambda)
% DAMPED_PINV simply computes the dampened pseudoinverse of the provided
%   matrix.
%   A general approach to solving the problem of discontinuity of the 
%   pseudoinverse solution at a singular configuration is to use the damped
%   least-squares method, which is known as the Levenberg-Marquardt 
%   stabilization method. This solution minimizes 
%   ||Xdot - J*qdot||^2 + lambda^2*||qdot||^2
%   with lambda constant.
%
%   INPUT
%   mat    - matrix to invert
%   lambda - dampening factor

    n = size(mat,1);
    core = mat * mat' + lambda^2 * eye(n);
    dpinv = mat' * inv(core);
end
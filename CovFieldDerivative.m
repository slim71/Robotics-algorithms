function Lfomegax = CovFieldDerivative(f, omega, vars)
% COVFIELDDERIVATIVE represents the directional derivative of the covector field
%   omega along the vector field f.
%   Notation: L_f omega(x) = f(x)^T *(domega^T/dx)^T + omega(x) * df/dx

    assert(iscolumn(f) && ~isscalar(f), "f must be a column vector " + ...
        "(it's a vector field)");
    assert(isrow(omega) && ~isscalar(omega), "g must be a row vector " + ...
        "(it's a covector field)");

    Lfomegax = f' * jacobian(omega', vars)' + omega * jacobian(f, vars);
end

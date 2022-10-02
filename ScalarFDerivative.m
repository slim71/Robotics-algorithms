function Lflambdax = ScalarFDerivative(vfield, lambda, vars)
% LIEBRACKET computes .
%   This represents the directional derivative of the scalar function lambda 
%   along the vector field f.

    assert(iscolumn(vfield) && ~isscalar(vfield), "vfield must be a column vector " + ...
        "(it's a vector field)");
    assert(isscalar(lambda), "lambda must be a scalar function");

    Lflambdax = jacobian(lambda, vars) * vfield;
end

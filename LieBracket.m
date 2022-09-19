function Lfgx = LieBracket(f, g, vars)
% LIEBRACKET computes the Lie-Bracket between vector fields f and g w.r.t. vars.
%   Notation: L_f g(x) = [f(x), g(x)] = ad_f g(x) = dg/dx * f(x) - df/dx * g(x)
%   This represents the directional derivative of the vector field g along
%   the vector field f.

    assert(iscolumn(f) && ~isscalar(f), "f must be a column vector " + ...
        "(it's a vector field)");
    assert(iscolumn(g) && ~isscalar(g), "g must be a column vector " + ...
        "(it's a vector field)");

    Lfgx = jacobian(g, vars) * f - jacobian(f, vars) * g;
end

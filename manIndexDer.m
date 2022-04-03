function varargout = manIndexDer(qq, jacobian, q_sym) 
    dj = jacobian(q_sym);
    for i=[1:size(q_sym,1)]
        der_J(:,:,i) = diff(dj,q_sym(i));
    end

    a = matlabFunction(der_J(:,:,1));
    a(qq(1), qq(2), qq(3), qq)
end
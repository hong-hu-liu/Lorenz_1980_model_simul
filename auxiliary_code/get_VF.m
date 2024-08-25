function f = get_VF(u_cur, G, Idx1, Idx2, A, Forcing)

% compute the vector field at a given state u_cur,
% where the vector field is of the form: A*u_cur + B(u_cur,u_cur) + Forcing,
% with B being a quadratic form defined by G. 

    nonlin = zeros(9,1);
    % compute the nonlinearity
    for k = 1:9
        temp = u_cur(Idx1{k}).*u_cur(Idx2{k});
        nonlin(k) = temp'*G{k};
    end   
    
    f = A*u_cur + nonlin + Forcing;
    
    
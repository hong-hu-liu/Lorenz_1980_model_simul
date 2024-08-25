% integrate the 9D Lorenz model using RK4  

function [u,PAR_Lorenz9D] = int_Lorenz9D(alpha, F1, u0, Ntmax, dt)

PAR_Lorenz9D = get_par_Lorenz9D(alpha, F1);


u = zeros(9,Ntmax);
u(:,1) = u0;

G = PAR_Lorenz9D.G; % coefficients for the quadratic terms 
Idx1 = PAR_Lorenz9D.Idx1;
Idx2 = PAR_Lorenz9D.Idx2;
A = PAR_Lorenz9D.A; % linear part
Forcing = PAR_Lorenz9D.Forcing;


fprintf('Solving the Lorenz-80 model using RK4...\n');
for i = 2:Ntmax

    if mod(i,5e4) == 0
        fprintf('i = %d out of %d...\n',i,Ntmax);
    end

    u_cur = u(:,i-1);

    %%% RK4
    k1 = dt*get_VF(u_cur, G, Idx1, Idx2, A, Forcing);
    k2 = dt*get_VF(u_cur+ 0.5*k1, G, Idx1, Idx2, A, Forcing);
    k3 = dt*get_VF(u_cur+ 0.5*k2, G, Idx1, Idx2, A, Forcing);
    k4 = dt*get_VF(u_cur+ k3, G, Idx1, Idx2, A, Forcing);
    u(:,i) = u_cur + (k1 + 2*k2 + 2*k3 + k4)./6;

    if isnan(u(1,i))
        warning(['Solution blows up at step ',num2str(i)]);
        return;
    end
end
fprintf('End of solving the Lorenz-80 model.\n\n');

    
    
% This code simulates the 9-dimensional ODE system introduced by Lorenz in
% 1980, as a Fourier spectral truncation of the Primitive Equations

% The model equations are those given by Eqns (33)-(35) in
% [L80] E. N. Lorenz (1980): Attractor sets and quasi-geostrophic equilibrium. J Atmos Sci 37, 1685–1699.
% We refer to this model as the L80 model.

% It is solved by using a fourth-order Runge–Kutta method (RK4) scheme. 

% See also Eqns. (2.1a)-(2.1c) in:
% [CLM17] M. D. Chekroun, H. Liu, and J. C. McWilliams (2017): The emergence of fast oscillations in a reduced 
% primitive equation model and its implications for closure theories. Comput. Fluids 151, 3–22.
% https://doi.org/10.1016/j.compfluid.2016.07.005


close all;
clear;
addpath('./auxiliary_code');

%--------------------------------

% Parameter regimes:

% The regimes provided below are the two regimes analyzed in:

% [CLSM24] M. D. Chekroun, H. Liu, K. Srinivasan, and J. C. McWilliams (2024): The high-frequency and rare events 
% barriers to neural closures of atmospheric dynamics. J. Physics Complexity, 5, 025004.
% https://iopscience.iop.org/article/10.1088/2632-072X/ad3e59

% and in

% [CLM21] M.D. Chekroun, H. Liu, and J. McWilliams (2021): Stochastic rectification of fast oscillations on slow manifold closures.
% Proc. Natl. Acad. Sci. (PNAS), 118(48):e2113650118. https://doi.org/10.1073/pnas.2113650118

% The system parameters are given by Eq. (2.3) in [CLM17], except for  the constant forcing F1.

% The latter parameter is a control parameter. Different values of F1 lead to different dynamics ranging from periodic/quasi-periodic solutions to
% slow Lorenz63-like chaos and to chaos with explosive fast oscillations as analyzed in [CLM17] (see Figure 6 therein).
% Note that the change of variables given by [CLM17, Eq. (2.4)] is meant to identify a parameter espilon related to the Rossby number.
% Changing epsilon corresponds also to changing F_1. The model simulated below is written in
% the original formulation of Lorenz in 1980.


par_id = 1; % High-low frequency (HLF) regime analyzed in [CLSM24] 
%par_id = 2; % Slow chaos regime analyzed in [CLSM24]

if par_id == 1
    F1 = 0.3027; 
elseif par_id == 2
    F1 =  0.0697;
else
    error('Need to specify value of F1.'); 
end

alpha = 3;  % This is the parameter a3 in the L80 model. 

%-----------------------
% set up time step size and total number of time steps: 
freq  =1;
dtFactor = 1/(2*freq)*1E-1; % A factor used to adjust the time step size relative to dt=1/12 (see next line) 
dt = 1/12*dtFactor; % dt = 1/12 corresponds to 15 minutes according to Lorenz's original L80 paper. 

Ntmax = 5e5; %  total number of time steps 

T = (Ntmax-1)/(4*24)*dtFactor; % max time of integration in days.  
tt = 0:0.25/24*dtFactor:T; % time vector corresponds to the time discretization used (unit: day)
%-----------------------


%---initial data----------- 
%The initial data is chosen to be close to the Hadley fixed point; see [CLM17, Sec. 2.3] and  [CLM17, Eq. (2.3)]

u0 = zeros(9,1);  % The first three components of u correspond to x1, x2, x3, respectively, 
                           % the components 4-6 of u correspond to y1, y2, y3, respectively, and
                           % the components 7-9 of u correspond to z1, z2, z3, respectively. 

a1 = 1;     
nu0 = 1/48;
g0 = 8;
u0(4) = F1/(a1*nu0*(1+a1*g0));  % initial data for y1
u0(1) = -nu0*a1*u0(4);  % initial data for x1
u0(7) = u0(4);  % initial data for z1
u0(5) = -10^-5;  % initial data for y2
u0(8) = 10^-5;   % initial data for z2
%-------------------------------


[u,PAR_Lorenz9D] = int_Lorenz9D(alpha, F1, u0, Ntmax, dt);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% plot the embedded attractors
is_attr = round(Ntmax/2); 
ie_attr = Ntmax;
ie_plot=ie_attr;
y = u(4:6,:);
freq_attr = 5;

if par_id == 1
fm=1.5;
fM=0.9;
is_plot=floor(fm*is_attr);
ie_plot=floor(fM*ie_attr);
    figure('position', [20 300 1750 600]);
    subplot(121)
    plot(tt(is_plot:ie_plot), y(2,is_plot:ie_plot),'k-','LineWidth',1.6)
    grid on
    set(gca,'fontsize',18,'fontweight','b','LineWidth',1.2)
     ylabel('$y_2(t)$','interpreter','latex','fontsize',26);
     xlabel('$t$ (days)','interpreter','latex','fontsize',26);
      title('Time series','Fontsize',20,'FontWeight','Bold');
    xlim([tt(is_plot) tt(ie_plot)])
    ym=1.05*min(y(2,is_plot:ie_plot));
    yM=1.05*max(y(2,is_plot:ie_plot));
    ylim([ym yM])
    %------------------------------------------------------------%
    subplot(122)
    plot(y(3,is_attr:freq_attr:ie_attr), y(2,is_attr:freq_attr:ie_attr),'k','LineWidth',1.2)
     set(gca,'fontsize',18,'fontweight','b','LineWidth',1.5)
    xlabel('$y_3$','interpreter','latex','fontsize',26);
    ylabel('$y_2$','interpreter','latex','fontsize',26);
    title('Projected Attractor','Fontsize',22,'FontWeight','Bold');
    grid on;
else
    
fm=1;
is_plot=floor(fm*is_attr);
    figure('position', [20 300 1750 600]);
    subplot(121)
    plot(tt(is_plot:ie_plot), y(3,is_plot:ie_plot),'k-','LineWidth',2)
    grid on
     set(gca,'fontsize',18,'fontweight','b','LineWidth',1.2)
      ylabel('$y_3(t)$','interpreter','latex','fontsize',26);
     xlabel('$t$ (days)','interpreter','latex','fontsize',26);
           title('Time series','Fontsize',22,'FontWeight','Bold');
    xlim([tt(is_plot) tt(ie_plot)])
    ym=1.05*min(y(3,is_plot:ie_plot));
    yM=1.05*max(y(3,is_plot:ie_plot));
    ylim([ym yM])
    %------------------------------------------------------------%
    subplot(122)
    plot(y(3,is_attr:freq_attr:ie_attr), y(1,is_attr:freq_attr:ie_attr),'k-','LineWidth',2)
    set(gca,'fontsize',18,'fontweight','b','LineWidth',1.2)
    xlabel('$y_3$','interpreter','latex','fontsize',26);
    ylabel('$y_2$','interpreter','latex','fontsize',26);
    title('Projected Attractor','Fontsize',20,'FontWeight','Bold');
    grid on;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

return; 



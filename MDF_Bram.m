%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Molecular Dynamics for Fluids
% Author: Bram ter Huurne
% Course: APIE
% Date: 31/01/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% c's
clear all;
close all;
clc

% plot progression?
plotting = true;

% initialize particles
N = 150;
Eps = 100;
mass = 1;
Sigma = 1;

% initialize wall
L = 10;
Eps_wall = 200;
Sigma_wall = 1;

% initialize time
N_steps = 1E4;
dt = 5e-4;
mvm_size = [10,0];

% initialize energies, temp, pressure
E_pot = zeros(N_steps,1)';
E_kin = zeros(N_steps,1)';
T = zeros(N_steps,1)';
P = zeros(N_steps,1)';
P_law = zeros(N_steps,1)';

% initialize all particles
x = (L-2*Sigma_wall) * (rand(N, 2) - 0.5);
x_old = x;
v = zeros(size(x));
v_abs_sqr = zeros([length(x),N_steps]); %v;
v_max = 100; % for Maxwell-Boltzmann
v_arr = 1:1:v_max;

% init
f_ij = zeros(N,N,2);
phi_ij = zeros(N,N);
kB = 1; % physconst('Boltzmann');

% create permutation landscape
[i,j] = meshgrid(1:N, 1:N);

%%

if plotting
    init_figure;
else
    % init waitbar
    f = waitbar(0,'Starting simulation...');
end

% start integration
for step = 1:N_steps
    
    % particle distances
    r_ij(:, :, 1) = reshape(x(i,1) - x(j,1), N, N);
    r_ij(:, :, 2) = reshape(x(i,2) - x(j,2), N, N);
    r_ij_abs = sqrt(r_ij(:,:,1).^2 + r_ij(:,:,2).^2);
    
    % forces on particles
    f_ij = calculate_force(r_ij,r_ij_abs,Sigma,Eps);
    f_i = squeeze(sum(f_ij .* ~eye(N),1));
    
    % wall interaction
    r_ij_wall(:, 1) = x(:,1) - L/2;
    r_ij_wall(:, 2) = zeros(size(r_ij_wall(:, 1)));
    r_ij_wall_abs(:,:,1) = abs(r_ij_wall(:, 1));
    f_i_wall = calculate_force(r_ij_wall,r_ij_wall_abs(:,:,1),Sigma_wall,Eps_wall);
    
    r_ij_wall(:, 1) = x(:,1) + L/2;
    r_ij_wall(:, 2) = zeros(size(r_ij_wall(:, 1)));
    r_ij_wall_abs(:,:,2) = abs(r_ij_wall(:, 1));
    f_i_wall = f_i_wall + calculate_force(r_ij_wall,r_ij_wall_abs(:,:,2),Sigma_wall,Eps_wall);
    
    r_ij_wall(:, 2) = x(:,2) - L/2;
    r_ij_wall(:, 1) = zeros(size(r_ij_wall(:, 2)));
    r_ij_wall_abs(:,:,3) = abs(r_ij_wall(:, 2));
    f_i_wall = f_i_wall + calculate_force(r_ij_wall,r_ij_wall_abs(:,:,3),Sigma_wall,Eps_wall);
    
    r_ij_wall(:, 2) = x(:,2) + L/2;
    r_ij_wall(:, 1) = zeros(size(r_ij_wall(:, 2)));
    r_ij_wall_abs(:,:,4) = abs(r_ij_wall(:, 2));
    f_i_wall = f_i_wall + calculate_force(r_ij_wall,r_ij_wall_abs(:,:,4),Sigma_wall,Eps_wall);
    
    % total force 
    f_i = f_i + f_i_wall;
       
    % calculate energies
    phi_ij = Eps .* ( Sigma - r_ij_abs).^2 .*(r_ij_abs <= Sigma) .* ~eye(N);
    phi_wall = Eps_wall .* ( Sigma_wall - r_ij_wall_abs).^2 .*(r_ij_wall_abs <= Sigma_wall);
    E_pot(step) = 0.5 * squeeze(sum(sum(phi_ij))) + squeeze(sum(sum(phi_wall)));
    v_abs_sqr(:,step) = v(:,1).^2 + v(:,2).^2;
    E_kin(step) = sum(0.5 .* mass .* v_abs_sqr(:,step));
    
    % temperature
    T(step) = E_kin(step) / (N * kB);
    
    % pressure: wall
    f_dot_n = sum(abs(f_i_wall),2);
    P_wall = 1/(4*L) * mean(f_dot_n);
    
    % pressure: virial
    f_ij_dot_r_ij = sum(r_ij.*f_ij,3);
    P_virial = N*kB*T(step)/L^2 + 1/(2*L^3) * mean(f_ij_dot_r_ij(:));
    P(step) = P_wall + P_virial;
    
    % total pressure
    P_law(step) = kB* N * T(step) / L^2;
    
    
    % integrate positions with Verlet
    tmp = x;
    [x,v] = verlet(x,x_old,f_i,dt,mass);
    x_old = tmp;
    
    avg_T = movmean(T,mvm_size);
    avg_P = movmean(P,mvm_size);
    avg_P_law = movmean(P_law,mvm_size);
    avg_v_abs_sqr = movmean(v_abs_sqr,mvm_size,2);
    
    % Maxwell-Boltzmann curve
    MWBM_curve = N .* v_arr .* exp(-v_arr.^2 .* mass./(2*kB.*avg_T(step))) ...
            .* mass ./ (kB .* T(step));

    if plotting && rem(step,round(N_steps/100)) == 0 % && step > N_steps/2
        N_in = sum(sqrt(x(:,1).^2 + x(:,2).^2)<(L/sqrt(2)));
        sgtitle(sprintf('Step: %d of %d  -  %d of %d particles in grid',step,N_steps,N_in,N));
        
        
        set(h1,'XData',x(:,1),'YData',x(:,2));
        xlim([-L/2,L/2]);
        ylim([-L/2,L/2]);
        
        set(h2_1,'YData',E_kin);
        set(h2_2,'YData',E_pot);
        set(h2_3,'YData',E_pot+E_kin);
        
        set(h3_1,'YData',avg_T);
        set(h4_1,'YData',avg_P);
        set(h4_2,'YData',avg_P_law);
        
        set(h5_1,'Data',avg_v_abs_sqr(:,step))
        set(h5_2,'YData',MWBM_curve)
        pause(0.01)
        
    elseif ~plotting && rem(step,round(N_steps/100)) == 0
        waitbar(step/N_steps,f,sprintf('Step: %d of %d',step,N_steps));
    end
end

if ~plotting
    delete(f)
    plot_energy(E_kin,E_pot)
end

g = radial_dist(r_ij_abs,1,L);

curdate = datestr(datetime,'yyyy_mm_dd_HH_MM_SS');
savefig(h,['results//',curdate,'.fig'])
save([curdate,'.mat'],'g','E_kin','E_pot','avg_T','avg_P','avg_P_law','avg_v_abs_sqr','MWBM_curve'



function f_ij = calculate_force(r_ij,r_ij_abs,Sigma,Eps)
    f_ij = 2.*Eps.*(Sigma - r_ij_abs).*(r_ij_abs <= Sigma) .* r_ij ./ r_ij_abs;
    f_ij(isnan(f_ij)) = 0;
end

function [x_out,v_out] = verlet(x_in,x_in_old,f_i,dt,mass)
    x_out = 2 .* x_in - x_in_old + f_i.*dt^2./mass;
    v_out = (x_out - x_in) ./ dt;
end

function g = radial_dist(r_ij_abs,dr,L)
    r_max = sqrt(2)*L;
    rr = 0:dr:r_max;
    g = zeros(size(rr));
    for idx = 1:length(rr)
        r = rr(idx);
        g(idx) = sum(sum(r_ij_abs(r_ij_abs >= r & r_ij_abs < r+dr)));
    end
    
    figure; 
    plot(rr,g);
    xlabel('r (m)')
    ylabel('g(r)')    
end

function [] = plot_energy(E_kin,E_pot)
    t_steps = 1:length(E_kin);
    
    figure;
    plot(t_steps,E_kin);
    hold on;
    plot(t_steps,E_pot);
    plot(t_steps,E_pot+E_kin);
    legend('E_{kin}','E_{pot}','E_{tot}');
end


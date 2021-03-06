%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Molecular Dynamics for Fluids
% Author: Bram ter Huurne
% ID: s1491784
% Course: APIE
% Date: 31/01/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% c's
clear all;
% close all;
clc

% plot progression?
plotting = true;

% particle properties
N_arr = [150];
% N_arr = [10,20,50,100,200,500,1000];

% init saving arrays
avg_T_arr = zeros(size(N_arr));
avg_P_arr = zeros(size(N_arr));
avg_P_law_arr = zeros(size(N_arr));

for idx = 1:length(N_arr)
    N = N_arr(idx);
    Eps = 100;
    mass = 1;
    Sigma = 1;

    % wall properties
    L = 10;
    Eps_wall = 500;
    Sigma_wall = 1;

    % initialize time
    N_steps = 1E4; %2E3;
    dt = 1e-4; 
    mvm_size = [100,0];

    % initialize energies, temp, pressure
    E_pot = zeros(N_steps,1)';
    E_kin = zeros(N_steps,1)';
    Momentum_tot = zeros(N_steps,1)';
    T = zeros(N_steps,1)';
    P = zeros(N_steps,1)';
    P_law = zeros(N_steps,1)';

    % initialize all particles
    x = zeros([N,2,N_steps]);
    x(:,:,1) = (L-2*Sigma_wall) * (rand(N, 2) - 0.5);
    v = zeros(size(x(:,:,1)));
    v_abs_sqr = zeros([N,N_steps]);
    v_max = 40; % static to plot Maxwell-Boltzmann
    v_arr = 1:1:v_max;

    % init
    f_ij = zeros(N,N,2);
    phi_ij = zeros(N,N);
    r_ij = zeros(N,N,2);
    r_ij_wall = zeros(N,2);
    r_ij_wall_abs = zeros(N,4);
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
    for step = 1:N_steps-1

        % particle distances
        r_ij(:, :, 1) = reshape(x(i,1,step) - x(j,1,step), N, N);
        r_ij(:, :, 2) = reshape(x(i,2,step) - x(j,2,step), N, N);
        r_ij_abs = sqrt(r_ij(:,:,1).^2 + r_ij(:,:,2).^2);

        % forces on particles
        f_ij = calculate_force(r_ij,r_ij_abs,Sigma,Eps);
        f_i = squeeze(sum(f_ij .* ~eye(N),1));

        % wall interactions
        r_ij_wall(:, 1) = x(:,1,step) - L/2;
        r_ij_wall(:, 2) = zeros(size(r_ij_wall(:, 1)));
        r_ij_wall_abs(:,1) = abs(r_ij_wall(:, 1));
        f_i_wall = calculate_force(r_ij_wall,r_ij_wall_abs(:,1),Sigma_wall,Eps_wall);

        r_ij_wall(:, 1) = x(:,1,step) + L/2;
        r_ij_wall(:, 2) = zeros(size(r_ij_wall(:, 1)));
        r_ij_wall_abs(:,2) = abs(r_ij_wall(:, 1));
        f_i_wall = f_i_wall + calculate_force(r_ij_wall,r_ij_wall_abs(:,2),Sigma_wall,Eps_wall);

        r_ij_wall(:, 2) = x(:,2,step) - L/2;
        r_ij_wall(:, 1) = zeros(size(r_ij_wall(:, 2)));
        r_ij_wall_abs(:,3) = abs(r_ij_wall(:, 2));
        f_i_wall = f_i_wall + calculate_force(r_ij_wall,r_ij_wall_abs(:,3),Sigma_wall,Eps_wall);

        r_ij_wall(:, 2) = x(:,2,step) + L/2;
        r_ij_wall(:, 1) = zeros(size(r_ij_wall(:, 2)));
        r_ij_wall_abs(:,4) = abs(r_ij_wall(:, 2));
        f_i_wall = f_i_wall + calculate_force(r_ij_wall,r_ij_wall_abs(:,4),Sigma_wall,Eps_wall);

        % total force 
        f_i = f_i + f_i_wall;

        % calculate energies
        phi_ij = Eps .* ( Sigma - r_ij_abs).^2 .*(r_ij_abs <= Sigma) .* ~eye(N);
        phi_wall = Eps_wall .* ( Sigma_wall - r_ij_wall_abs).^2 .*(r_ij_wall_abs <= Sigma_wall);
        E_pot(step) = 0.5 * squeeze(sum(sum(phi_ij))) + squeeze(sum(sum(phi_wall)));
        v_abs_sqr(:,step) = sqrt(v(:,1).^2 + v(:,2).^2);
        E_kin(step) = sum(0.5 .* mass .* v_abs_sqr(:,step).^2);
        Momentum_tot(step) = sum(v_abs_sqr(:,step) * mass);

        % temperature
        T(step) = E_kin(step) / (N * kB);

        % pressure: wall
        f_dot_n = sum(abs(f_i_wall),2);
        P_wall = 1/(4*L) * mean(f_dot_n);

        % pressure: virial
        f_ij_dot_r_ij = sum(r_ij.*f_ij,3);
        P_virial = N*kB*T(step)/L^2 + 1/(2*L^3) * mean(f_ij_dot_r_ij(:));
        
        % total pressure
        P(step) = P_wall + P_virial;

        % total pressure by ideal gas law
        P_law(step) = kB* N * T(step) / L^2;


        % integrate positions with Verlet
        if step == 1
            [x(:,:,step+1),v] = verlet(x(:,:,1),x(:,:,1),f_i,dt,mass);
        else
            [x(:,:,step+1),v] = verlet(x(:,:,step),x(:,:,step-1),f_i,dt,mass);
        end

        % calculate averages
        avg_T = movmean(T,mvm_size);
        avg_P = movmean(P,mvm_size);
        avg_P_law = movmean(P_law,mvm_size);
        avg_v_abs_sqr = movmean(v_abs_sqr,mvm_size,2);

        % Maxwell-Boltzmann curve
        MWBM_curve = 2 * N .* v_arr .* exp(-v_arr.^2 .* mass./(2*kB.*avg_T(step))) ...
                .* mass ./ (kB .* T(step));

        if plotting && rem(step+1,round(N_steps/100)) == 0
            % show if particles have left domain
            N_in = sum(sqrt(x(:,1,step+1).^2 + x(:,2,step+1).^2)<(L/sqrt(2)));
            sgtitle(sprintf('Step: %d of %d  -  %d of %d particles in grid',step+1,N_steps,N_in,N));

            set(h1,'XData',x(:,1,step+1),'YData',x(:,2,step+1));
            xlim([-L/2,L/2]);
            ylim([-L/2,L/2]);

            set(h2_1,'YData',E_kin);
            set(h2_2,'YData',E_pot);
            set(h2_3,'YData',E_pot+E_kin);

            set(h6_1,'YData',Momentum_tot);

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
    % save figures
    curdate = datestr(datetime,'yyyy_mm_dd_HH_MM_SS');
    scaling=4; set(gcf, 'Color', 'none');
    export_fig (['results//',curdate,'.png'],['-m ',num2str(scaling)])
    savefig(h,['results//',curdate,'.fig']);
    
    % plot radial distribution function and self-diffusion coefficient
    gr = radial_dist(r_ij_abs,1,L,N);
    [gt,D] = mean_square_displacement(x,N_steps,dt);
    
    % save results
    save(['results//',curdate,'.mat'],'gr','gt','D','E_kin','E_pot','avg_T','avg_P','avg_P_law','avg_v_abs_sqr','MWBM_curve');

    % save averages to array
    avg_T_arr(idx) = avg_T(end-1);
    avg_P_arr(idx) = avg_P(end-1);
    avg_P_law_arr(idx) = avg_P_law(end-1);
end

% disp(N_arr);
% disp(avg_T_arr);
% disp(avg_P_arr);
% disp(avg_P_law_arr);

function f_ij = calculate_force(r_ij,r_ij_abs,Sigma,Eps)
    f_ij = 2.*Eps.*(Sigma - r_ij_abs).*(r_ij_abs <= Sigma) .* r_ij ./ r_ij_abs;
    f_ij(isnan(f_ij)) = 0;
end

function [x_out,v_out] = verlet(x_in,x_in_old,f_i,dt,mass)
    x_out = 2 .* x_in - x_in_old + f_i.*dt^2./mass;
    v_out = (x_out - x_in) ./ dt;
end

function g = radial_dist(r_ij_abs,dr,L,N)
    r_max = sqrt(2)*L;
    rr = 0:dr:r_max;
    g = zeros(size(rr));
    for idx = 1:length(rr)
        r = rr(idx);
        expected = 2*pi*N*r*dr/L^2;
        g(idx) = sum(sum(r_ij_abs(r_ij_abs >= r & r_ij_abs < r+dr)))/expected;
    end
    
    figure; 
    plot(rr,g);
    xlabel('r (m)')
    ylabel('g(r)')    
end


function [g,Ds] = mean_square_displacement(x,N_steps,dt)
    tau = 1000;
    tt = tau:N_steps;
    g = zeros(size(tt));
    for idx = 1:length(tt)-1
        g(idx) = mean((x(:,1,tau+idx)-x(:,1,tau)).^2 + ...
            (x(:,2,tau+idx)-x(:,2,tau)).^2);
    end
    
    D = g./(2.*2.*tt.*dt);
    Ds = D(end);
    
    figure; 
    subplot(1,2,1)
    plot(tt,g);
    hold on;
    plot(tt,2.*2.*tt.*dt*0.01);
    plot(tt,2.*2.*tt.*dt*0.02);
    plot(tt,2.*2.*tt.*dt*0.1);
    plot(tt,2.*2.*tt.*dt*0.2);
    plot(tt,2.*2.*tt.*dt*1);
    xlabel('t (s)')
    legend('g(t)','D_s=0.01','D_s=0.02','D_s=0.1','D_s=0.2','D_s=1')  
    hold off;
    subplot(1,2,2)
    plot(tt,D);
    xlabel('t (s)')
    ylabel('D_s')  
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


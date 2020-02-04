% init figure
h = figure('WindowState','maximized');%'units','normalized','outerposition',[0 0 1 1])

t_steps = (1:N_steps).*dt;
h_sp2 = subplot(3,3,3);
h2_1 = plot(t_steps,E_kin);
hold on;
h2_2 = plot(t_steps,E_pot);
h2_3 = plot(t_steps,E_pot+E_kin);
legend('E_{kin}','E_{pot}','E_{tot}');
xlabel('Time (s)')
ylabel('Energy (J)')
hold off;

h_sp3 = subplot(3,3,6);
h3_1 = plot(t_steps,T);
xlabel('Time (s)')
ylabel('Temperature (K)')

h_sp4 = subplot(3,3,9);
h4_1 = plot(t_steps,P,'+-','MarkerSize',3);
hold on;
h4_2 = plot(t_steps,P_law,'-','MarkerSize',3);
legend('P','P_{ideal}')
xlabel('Time (s)')
ylabel('Pressure (Pa)')

h_sp5 = subplot(3,3,[7,8]);
h5_1 = histogram(v_abs_sqr,v_max,'BinLimits',[0,v_max]);
hold on;
h5_2 = plot(v_arr,zeros(size(v_arr)));
legend('Maxwell','Maxwell_{ideal}')
xlabel('Velocity (m/s)')
ylabel('Count')


h_sp1 = subplot(3,3,[1,2,4,5]);
h1 = scatter(x(:,1),x(:,2),'Ob') ;
hold on;
box on;
sgtitle(sprintf('Step: %d of %d',0,N_steps));
xlim([-L/2,L/2]);
ylim([-L/2,L/2]);
axis equal;
hold off;
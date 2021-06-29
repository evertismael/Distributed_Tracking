function show_measurements(fig,t_vect, deltas_mean, deltas_var)
figure(fig)
    subplot(2,1,1);
    plot(t_vect,deltas_mean)
    title('Range Mean'); legend(string(1:12));
    xlabel('t [s]'); ylabel('\delta [m]'); grid on;

    subplot(2,1,2);
    semilogy(t_vect, sqrt(deltas_var))
    ylim([10^-1, 10^3]);
    title('Range Standart Deviation');
    xlabel('t [s]'); ylabel('\sigma [m]'); grid on;
end
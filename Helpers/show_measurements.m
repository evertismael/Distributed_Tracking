function show_measurements(fig,t_vect, deltas_mean, deltas_var)
figure(fig)
    subplot(2,1,1);
    plot(t_vect,deltas_mean)
    title('\delta'); legend();

    subplot(2,1,2);
    semilogy(t_vect,deltas_var)
    title('\Sigma'); legend();
end
function show_target(fig, target)
figure(fig)
    subplot(2,3,1);
    plot(target.t_vect,target.history(1,:))
    title('x')

    subplot(2,3,4);
    plot(target.t_vect,target.history(2,:))
    title('vx')

    subplot(2,3,2);
    plot(target.t_vect,target.history(3,:))
    title('y')

    subplot(2,3,5);
    plot(target.t_vect,target.history(4,:))
    title('vy')

    subplot(2,3,6);
    hold on;
    plot(target.history(1,:),target.history(3,:),'.b');
    plot(target.history(1,1),target.history(3,1),'or');
    
    title('xy-plane')
    xlim([0,60]); ylim([0,60]);
end
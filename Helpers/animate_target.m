function animate_target(fig, target, xy_dim)
    figure(fig)
    title('xy-plane')
    axis(xy_dim); grid on;
    hold on;
    scene = Params.get_scene();
    scatter(scene.bx(1,:),scene.bx(2,:),'^');
    circle(45,45,5); % roundabout:
    road_borders();
    for idx = 1:size(target.history(3,:),2)-1
        plot(target.history(1,idx),target.history(3,idx),'.b');
        pause(0.1);
    end  
end


function h = circle(x,y,r)
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h = plot(xunit, yunit);
end

function road_borders()
plot([55 55],[0 30],'k');
plot([35 35],[0 30],'k');

plot([55 55],[60 90],'k');
plot([35 35],[60 90],'k');

plot([0 30],[35 35],'k');
plot([0 30],[55 55],'k');

plot([60 90],[35 35],'k');
plot([60 90],[55 55],'k');

plot([30 35],[55 60],'k');
plot([30 35],[35 30],'k');
plot([55 60],[30 35],'k');
plot([55 60],[60 55],'k');

plot([45 45],[0 30],'k');
plot([45 45],[60 90],'k');

plot([0 30],[45 45],'k');
plot([60 90],[45 45],'k');
end
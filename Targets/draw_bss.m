function draw_bss()
    scene = Params.get_scene();
    scatter(scene.bx(1,:),scene.bx(2,:),'^');
    str = string(1:12);
    textscatter(scene.bx(1,:)+1,scene.bx(2,:)+1,str);
end
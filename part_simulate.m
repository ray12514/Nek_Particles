clear all
clc
close all
clf

scale=1;
for ii=001:3
    iname=sprintf('part%05d.3D',ii);
    iname2=sprintf('vel%05d.3D',ii);
    part = load(iname);
    vel = load(iname);
    h=plot3(part(:,1),part(:,2),part(:,3),'k.','markersize',10);
    hold on
    g=quiver3(part(:,1),part(:,2),part(:,3),...
              vel(:,1),vel(:,2),vel(:,3),scale);
     axis equal;axis ([-0.5 0.5 0.0 2.0 -0.5 0.5]);view([90,-90]);
     xlabel('x')
     ylabel('y')
     zlabel('z')
     drawnow
    hold off
% pause() 
%   F = getframe(fig);
%   aviobj = addframe(aviobj,F);
end
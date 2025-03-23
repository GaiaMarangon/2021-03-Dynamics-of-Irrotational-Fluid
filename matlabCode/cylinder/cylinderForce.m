clear all, close all


%-- PARAMETERS -------------------------------------------------------
% general
    U = -20;
    a = 2.0;
    %gamma = + 4*pi*a*U*0.5;   %positive for clockwise 
    %gamma = - 4*pi*a*U*1.0;   %negative for anticlockwise
    gamma = 0;                %zero for no circulation
%-- PARAMETERS -------------------------------------------------------


%-- SETTING DOMAIN --------------------------------------------------
 % building cylinder in the middle
    npoints = 80;
    interi=linspace(0,npoints,npoints)';
    theta=2*pi*interi/npoints;
    xtheta = a*cos(theta);
    ytheta = a*sin(theta);
%-- SETTING DOMAIN --------------------------------------------------


%-- FORCE FIELD ---------------------------------------------------
% pressure at surf
    p = ( 1 - (2*sin(theta)+gamma/(2*pi*a*U)).^2 )/2;
% force 
    fnorm = - p * 2*pi/(npoints-1) * a;
    fx = fnorm.*cos(theta);
    fy = fnorm.*sin(theta);
% plot
    figure(1)
    hold on 
    fill(xtheta,ytheta,'k','HandleVisibility','off')
    quiver(xtheta,ytheta,fx,fy,0,'color',[0.3010 0.7450 0.9330],'LineWidth',2); %light blue
    legend('Force Field');
    axis equal;
    xlabel('x');
    ylabel('y');
    xlim([-5,5]);
    ylim([-5,5]);
% title('Force Field');
    hold off    
%-- FORCE FIELD ---------------------------------------------------




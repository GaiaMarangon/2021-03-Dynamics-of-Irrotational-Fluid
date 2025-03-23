clear all, close all


%-- PARAMETERS --------------------------------------------------
    b = 5.0;
    e = 0.1;        % x shifting
    beta = pi/18;   % y shifting
    a = b*(e+1);
    radius = a/cos(beta);
    ycenter = a*tan(beta);
    xcenter = b*e;
%-- PARAMETERS --------------------------------------------------


%-- CYLINDER ----------------------------------------------------
% building points of cylinder
    npoints = 500;
    interi = linspace(0,npoints-1,npoints);
    theta = 2*pi*interi/(npoints-1);
    etheta = radius*cos(theta) + xcenter; 
    ntheta = radius*sin(theta) + ycenter; 
% plot cylinder
    figure(1)
    hold on
    axis equal
    % draw axis
        xdom=[-a-2,a+2];
        ydom=[-a-2,a+2];
        xaxx=linspace(xdom(1),xdom(2));
        yaxx=linspace(0,0);
        plot(xaxx,yaxx,'k-');
        xax=linspace(0,0);
        yax=linspace(ydom(1),ydom(2));
        plot(xax,yax,'k-');
    plot(etheta,ntheta,'b','LineWidth',2)
    plot(xcenter,ycenter,'xb','LineWidth',2)
    xlabel('\xi');
    ylabel('\eta');
    xlim(xdom);
    ylim(ydom);
    hold off
%-- CYLINDER ----------------------------------------------------


%-- AIRFOIL -----------------------------------------------------
% punti airfoil
    xair = etheta + (b^2*etheta)./(etheta.^2 + ntheta.^2);
    yair = ntheta - (b^2*ntheta)./(etheta.^2 + ntheta.^2);
% % punti chord
%     [minchord,iminchord] = min(xair);
%     [maxchord,imaxchord] = max(xair);
%     nchord = 100;
%     offchord = 1;
%     xchord = linspace(minchord,maxchord,nchord);
%     ychord = ones(1,nchord)*(min(yair)-offchord);
%     x1vertichord = ones(1,nchord/10)*maxchord;
%     x2vertichord = ones(1,nchord/10)*minchord;
%     y1vertichord = linspace(min(yair)-offchord,yair(iminchord),nchord/10);
%     y2vertichord = linspace(min(yair)-offchord,yair(imaxchord),nchord/10);
% plot airfoil 
    figure(2)
    hold on
    axis equal
    % draw axis
        xdom=[-2.5*b,2.5*b];
        ydom=[-2.5*b,2.5*b];
        xaxx=linspace(xdom(1),xdom(2));
        yaxx=linspace(0,0);
        plot(xaxx,yaxx,'k-');
        xax=linspace(0,0);
        yax=linspace(xdom(1),xdom(2));
        plot(xax,yax,'k-');
    xlabel('x');
    ylabel('y');
    xlim(xdom);
    ylim(ydom);
    % airfoil
        plot(xair,yair,'r','LineWidth',2)
%     % chord
%         plot(xchord,ychord,'b')
%         plot(x1vertichord,y1vertichord,':b');
%         plot(x2vertichord,y2vertichord,':b');
    hold off
%-- AIRFOIL -----------------------------------------------------

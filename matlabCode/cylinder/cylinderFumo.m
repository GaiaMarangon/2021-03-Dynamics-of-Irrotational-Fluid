clear all, close all

%-- PARAMETERS -------------------------------------------------------
% general
    U = -0.5;
    a = 2.0;
    %gamma = + 4*pi*a*U*0.25;   %positive for clockwise 
    %gamma = - 4*pi*a*U*1.0;   %negative for anticlockwise
    gamma = 0;                %zero for no circulation
% for streaklines
    nlines = 31;
    dt = 0.1;
    ntime = 1000;
    nsplit = 6;
    interv = floor(ntime/nsplit);
    off =floor(interv/20); % tra 0 e interv
%-- PARAMETERS -------------------------------------------------------


%-- SETTING DOMAIN --------------------------------------------------
% domain limits
    xdom = [-10,10];
    ydom = [-10,10];
 % building cylinder in the middle
    npoints = 100;
    interi=linspace(0,npoints,npoints)';
    theta=2*pi*interi/npoints;
    xtheta = a*cos(theta);
    ytheta = a*sin(theta);
%-- SETTING DOMAIN --------------------------------------------------



%-- STREAKLINES -----------------------------------------------------
% setting starting points
    dy = (ydom(2)-ydom(1))/nlines;
    y0 = [linspace(-dy*(nlines-1)/2,0,(nlines+1)/2)'; 
          linspace(dy,dy*(nlines-1)/2,(nlines-1)/2)']; 
    x0 = ones(nlines,1)*xdom(2);
% compute their time evolution
    x = [x0,zeros(nlines,ntime-1)];
    y = [y0,zeros(nlines,ntime-1)];
    for k=2:ntime
        for i=1:nlines
            r0 = sqrt( x(i,k-1)^2+y(i,k-1)^2 );
            o0 = atan2( y(i,k-1), x(i,k-1) );
            ur = +   U*(1-a^2/r0^2)*cos(o0);
            uo = - ( U*(1+a^2/r0^2)*sin(o0) + gamma/(2*pi*r0) );
            ux = ur*cos(o0) - uo*sin(o0);
            uy = ur*sin(o0) + uo*cos(o0);
            x(i,k) = x(i,k-1) + ux*dt;
            y(i,k) = y(i,k-1) + uy*dt;
            if (x(i,k)^2+ y(i,k)^2>a^2+0.01)
            else
                x(i,k)=-a-0.01;
            end
        end
    end
% plot
    % set figure
    figure(1)
    hold on
    axis equal
    xlim(xdom)
    ylim(ydom)
    fill(xtheta,ytheta,'k','HandleVisibility','off')
    % draw streaklines with proper color splitting
    for i=1:nlines
        if (mod( abs(i-(nlines+1)/2),5 ) == 0)
                plot(x(i,:),y(i,:),'-k','LineWidth',2)
        else
            for j=1:nsplit
                xplot = x(i,max([1,1+interv*(j-1)-off]):interv*j-off);
                yplot = y(i,max([1,1+interv*(j-1)-off]):interv*j-off);
                if (mod(j,2) == 0)
                    plot(xplot,yplot,'-b','LineWidth',2)
                else
                    plot(xplot,yplot,'-r','LineWidth',2)
                end
            end
        end
    end
%-- STREAKLINES -----------------------------------------------------

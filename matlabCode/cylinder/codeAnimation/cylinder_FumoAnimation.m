clear all, close all

% FRAMES SETTING
% Define number of frames
nr_fr = 80;
% Initialize matrix using 'moviein'
frames = moviein(nr_fr); 
% Save in avi file
v = VideoWriter('cylinder_Fumo1Animation.avi','Uncompressed AVI');
v.FrameRate = 8;
open(v)

%-- PARAMETERS -------------------------------------------------------
% general
    U = -0.5;
    a = 2.0;
    gamma = + 4*pi*a*U*0.25;   
% for streaklines
    nlines = 31;
    dt = 0.1;
    ntime = 1600;
    nsplit = 12;
    interv = floor(ntime/nsplit);
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
%-- STREAKLINES -----------------------------------------------------



%-- ITERATIONS FOR FRAMES - ANDATA -----------------------------------------------------------------
for k = 1 : nr_fr
    
    % CYLINDER FLOW 
            off =floor( (nr_fr-k)*2*interv/(nr_fr-1) ); % tra 0 e interv
        % set figure
            figure(1)
            clf
            hold on
            axis equal
            xlim(xdom)
            ylim(ydom)
            fill(xtheta,ytheta,'k','HandleVisibility','off')
        % draw streaklines with proper color splitting
            for i=1:nlines
                if (mod( abs(i-(nlines+1)/2),5 ) == 0)
                        plot(x(i,:),y(i,:),'-k','LineWidth',1)
                else
                    for j=1:nsplit
                        xplot = x(i,max([1,1+interv*(j-1)-off]):interv*j-off);
                        yplot = y(i,max([1,1+interv*(j-1)-off]):interv*j-off);
                        if (mod(j,2) == 0)
                            plot(xplot,yplot,'-b','LineWidth',1)
                        else
                            plot(xplot,yplot,'-r','LineWidth',1)
                        end
                    end
                end
            end
    % END CYLINDER FLOW
    frame = getframe(gcf);
    writeVideo(v,frame);
    if (k == nr_fr || k==1)     %first frame: keep same figure for some frames
        for h=1:10
            frame = getframe(gcf);
            writeVideo(v,frame);
        end
    end
end
%-- ITERATIONS FOR FRAMES - ANDATA -----------------------------------------------------------------

close(v)

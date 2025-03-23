clear all, close all

% FRAMES SETTING
% Define number of frames
nr_fr = 80;
% Initialize matrix using 'moviein'
frames = moviein(nr_fr); 
% Save in avi file
v = VideoWriter('airfoil_FumoP18Animation.avi','Uncompressed AVI');
v.FrameRate = 8;
open(v)

%-- PARAMETERS -------------------------------------------------------
% for shifting, rotating
    alpha = pi/18;  % rotating
    beta = pi/18;   % y shifting
    e = 0.1;        % xshifting
    toll = 0.02;
% general
    U = -0.5;
    b = 2.0;
    a = b*(e+1);
    radius = a/cos(beta);
    gamma = 4*pi*radius*U*sin(alpha+beta);
    ncenter = a*tan(beta);
    ecenter = b*e;
% for streaklines
    nlines = 41;
    dt = 0.1;
    ntime = 1200;
    nsplit = 6;
    interv = floor(ntime/nsplit);
%-- PARAMETERS -------------------------------------------------------


%-- SETTING DOMAIN --------------------------------------------------
% domain limits
    xdom = [-10,10];
    ydom = [-10,10];
    edom = [-15,15];
    ndom = [-15,15];
 % points of cylinder
    npoints = 100;
    interi=linspace(0,npoints,npoints)';
    theta=2*pi*interi/npoints;
    etheta = radius*cos(theta)+ecenter;
    ntheta = radius*sin(theta)+ncenter;
 % points of airfoil
    xxair = etheta + (b^2*etheta)./(etheta.^2 + ntheta.^2);
    yyair = ntheta - (b^2*ntheta)./(etheta.^2 + ntheta.^2);
    [xair,yair] = rotcart2cart(xxair,yyair,alpha);
%-- SETTING DOMAIN --------------------------------------------------



%-- STREAKLINES CYLINDER --------------------------------------------
% setting starting points
    dn = (ndom(2)-ndom(1))/nlines;
    n0 = [linspace(-dn*(nlines-1)/2,0,(nlines+1)/2)'; 
          linspace(dn,dn*(nlines-1)/2,(nlines-1)/2)']; 
    e0 = ones(nlines,1)*edom(2);
% compute their time evolution
    e = [e0,zeros(nlines,ntime-1)];
    n = [n0,zeros(nlines,ntime-1)];
    for k=2:ntime
        for i=1:nlines
            r0 = sqrt( (e(i,k-1)-ecenter)^2+(n(i,k-1)-ncenter)^2 );
            o0 = atan2( n(i,k-1)-ncenter, e(i,k-1)-ecenter );
            ur = +   U*(1-radius^2/r0^2)*cos(o0);
            uo = - ( U*(1+radius^2/r0^2)*sin(o0) + gamma/(2*pi*r0) );
            ux = ur*cos(o0) - uo*sin(o0);
            uy = ur*sin(o0) + uo*cos(o0);
            e(i,k) = e(i,k-1) + ux*dt;
            n(i,k) = n(i,k-1) + uy*dt;
            if ( (e(i,k)-ecenter)^2 + (n(i,k)-ncenter)^2>radius^2+0.01)
            else
                e(i,k)=-radius-0.01;
            end
        end
    end
%-- STREAKLINES CYLINDER --------------------------------------------


%-- STREAKLINES AIRFOIL --------------------------------------------
% rotate cartesian
    r = (e.^2 + n.^2 ).^(0.5);
    o = atan2(n,e) - alpha;
    ee = r.*cos(o);
    nn = r.*sin(o);
% apply Joukowsky to rotated cartesian
    xx = ee + (b^2*ee)./(ee.^2 + nn.^2);
    yy = nn - (b^2*nn)./(ee.^2 + nn.^2);
% rotate final view
    [x,y] = rotcart2cart(xx,yy,alpha);
%-- STREAKLINES AIRFOIL --------------------------------------------



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
            fill(xair,yair,'k','HandleVisibility','off')
            %fill(etheta,ntheta,'k','HandleVisibility','off')
            % draw streaklines with proper color splitting
            for i=1:nlines
                if (mod( abs(i-(nlines+1)/2),5 ) == 0)
                        plot(x(i,:),y(i,:),'-k','Linewidth',1)
                else
                    for j=1:nsplit
                        xplot = x(i,max([1,1+interv*(j-1)-off]):interv*j-off);
                        yplot = y(i,max([1,1+interv*(j-1)-off]):interv*j-off);
                        if (mod(j,2) == 0)
                            plot(xplot,yplot,'-b','Linewidth',1)
                        else
                            plot(xplot,yplot,'-r','Linewidth',1)
                        end
                    end
                end
            end
    % END CYLINDER FLOW
    
    frame = getframe(gcf);
    writeVideo(v,frame);
end
%-- ITERATIONS FOR FRAMES - ANDATA -----------------------------------------------------------------

close(v)



%----------------------------------------------------------------------    
% from adapted (rotated) cartesian to final cartesian
function [xgrid,ygrid] = rotcart2cart(xxgrid,yygrid,alpha)
    ss = ( xxgrid.^2 + yygrid.^2 ).^(1/2);
    tt = atan2(yygrid,xxgrid)+alpha;
    xgrid = ss.*cos(tt);
    ygrid = ss.*sin(tt);
end
%----------------------------------------------------------------------
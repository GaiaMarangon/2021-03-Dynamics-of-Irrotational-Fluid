clear all, close all


%-- PARAMETERS -------------------------------------------------------
% general
    U = -20.0;
    a = 2.0;
    %gamma = + 4*pi*a*U*0.5;   %positive for clockwise 
    %gamma = - 4*pi*a*U*1.0;   %negative for anticlockwise
    gamma = 0;                %zero for no circulation
% for domain
    dim = 101;
    dom = [-10,10];
    lessdim = 20;
%-- PARAMETERS -------------------------------------------------------


%-- SETTING DOMAIN --------------------------------------------------
% general cartesian grid
    xgrid = linspace(dom(1),dom(2),dim)';
    ygrid = linspace(dom(1),dom(2),dim)';
    [xgrid,ygrid] = meshgrid(xgrid,ygrid);
% general polar grid
    rgrid = sqrt( xgrid.^2 + ygrid.^2 );
    ogrid = atan2( ygrid,xgrid );
% coarser cartesian grid for velocity
    xcoarse = linspace(dom(1),dom(2),lessdim)';
    ycoarse = linspace(dom(1),dom(2),lessdim)';
    [xcoarse,ycoarse] = meshgrid(xcoarse,ycoarse);
% coarser polar grid for velocity
    rcoarse = sqrt( xcoarse.^2 + ycoarse.^2 );
    ocoarse = atan2( ycoarse,xcoarse );
% building cylinder in the middle
    npoints = 100;
    interi=linspace(0,npoints,npoints)';
    theta=2*pi*interi/npoints;
    xtheta = a*cos(theta);
    ytheta = a*sin(theta);
%-- SETTING DOMAIN --------------------------------------------------


%-- STREAMLINES AND EQUIPOTENTIAL LINES ------------------------------
% functions
    phi = U*( rgrid + a^2./rgrid ).*cos(ogrid) - gamma/(2*pi) * ogrid;
    psi = U*( rgrid - a^2./rgrid ).*sin(ogrid) + gamma/(2*pi) * log(rgrid/a);
% const values for contour
    maxphi = 600;   
    maxpsi = 600; 
    minphi = -600;
    minpsi = -600; 
    phivec = linspace(minphi,maxphi,40)';
    psivec = linspace(minpsi,maxpsi,41)';
% plot
    figure(1)
    hold on 
    fill(xtheta,ytheta,'k','HandleVisibility','off')
    contour(xgrid,ygrid,phi,phivec,'color',[0.0000 0.4470 0.7410],'LineWidth',1); %blue
    contour(xgrid,ygrid,psi,psivec,'color',[0.4660 0.6740 0.1880],'LineWidth',1); %green
    legend('Equipot. lines','Streamlines');
    axis equal;
    xlabel('x');
    ylabel('y');
    xlim(dom);
    ylim(dom);
    %title('Streamlines and Equipotential lines');
    hold off
%-- STREAMLINES AND EQUIPOTENTIAL LINES ------------------------------


%-- VELOCITY COMPONENTS ----------------------------------------------
% functions
    ucomp = zeros(lessdim,lessdim);
    vcomp = zeros(lessdim,lessdim);
    for i=1:lessdim
        for j=1:lessdim
            if ( rcoarse(i,j)<a )
                urcomp = 0;
                uocomp = 0;
            else
                urcomp = U*(1-a^2/rcoarse(i,j)^2)*cos(ocoarse(i,j));
                uocomp = -U*(1+a^2/rcoarse(i,j)^2)*sin(ocoarse(i,j)) - gamma/(2*pi*rcoarse(i,j));
            end
            ucomp(i,j) = urcomp*cos(ocoarse(i,j)) - uocomp*sin(ocoarse(i,j));
            vcomp(i,j) = urcomp*sin(ocoarse(i,j)) + uocomp*cos(ocoarse(i,j));
        end
    end
% plot
    figure(2)
    hold on 
    quiver(xcoarse,ycoarse,ucomp,vcomp,'color',[0.4940 0.1840 0.5560],'LineWidth',1); %violet
    fill(xtheta,ytheta,'k','HandleVisibility','off')
    legend('Velocity Field');
    axis equal;
    xlabel('x');
    ylabel('y');
    xlim(dom);
    ylim(dom);
    %title('Velocity Field');
    hold off
%-- VELOCITY COMPONENTS ----------------------------------------------


%-- PRESSURE FIELD ---------------------------------------------------
% function
    p = zeros(dim,dim);
    for i=1:dim
        for j=1:dim
                if (rgrid(i,j) < a-0.2)
                    p(i,j) = 0;
                else
                    p(i,j) = - ( a^4/rgrid(i,j)^4 - 2*a^2/rgrid(i,j)^2*cos(2*ogrid(i,j)) +...
                         ( gamma/(2*pi*rgrid(i,j)*U) )^2 + gamma/(pi*rgrid(i,j)*U)...
                         *(1+a^2/rgrid(i,j)^2)*sin(ogrid(i,j))  )/2;
                end
        end
    end
% const values for contour
    maxp = 0.5;
    midp = -0.9;
    minp = -7.5;
    pvec = [linspace(minp,midp,150)';
            linspace(midp,maxp,15)'];
% pressure plot
    figure(3)
    hold on 
    contourf(xgrid,ygrid,p,pvec,'LineColor','none');
    colormap(hsv)
    colorbar('westoutside');
    caxis([-2,maxp])
    fill(xtheta,ytheta,'k','HandleVisibility','off')
    legend('Pressure Field');
    axis equal;
    xlabel('x');
    ylabel('y');
    xlim(dom);
    ylim(dom);
    %title('Pressure Field');
    hold off
%-- PRESSURE FIELD ---------------------------------------------------

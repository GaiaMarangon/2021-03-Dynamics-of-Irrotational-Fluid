clear all, close all


%-- PARAMETERS --------------------------------------------------------
% general
    gamma = 20.0;   
% for domain
    dim = 101;
    dom = [-10,10];
    lessdim = 17;
%-- PARAMETERS --------------------------------------------------------


%-- SETTING DOMAIN ----------------------------------------------------
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
%-- SETTING DOMAIN ----------------------------------------------------


%-- STREAMLINES AND EQUIPOTENTIAL LINES -------------------------------
% functions
    phi = zeros(dim,dim);
    psi = zeros(dim,dim);
    for i=1:dim
        for j=1:dim
            if (rgrid(i,j) ~= 0)
                psi(i,j) = -gamma/(2*pi) * log(rgrid(i,j));
            else
                psi(i,j) = +gamma*1.0e14;
            end
            phi(i,j) = gamma/(2*pi) * ogrid(i,j);
        end
    end
% const values for contour
    maxphi = max(max(phi));     
    maxpsi = 8; 
    minphi = min(min(phi)); 
    minpsi = min(min(psi));
    phivec = linspace(minphi,maxphi,25)';
    psivec = linspace(minpsi,maxpsi,25)';
% plot
    figure(1)
    hold on 
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
%-- STREAMLINES AND EQUIPOTENTIAL LINES -------------------------------



%-- VELOCITY COMPONENTS -----------------------------------------------
% functions
    uocomp = gamma./(2*pi*rcoarse);
    urcomp = zeros(lessdim,lessdim);
    ucomp = urcomp.*cos(ocoarse) - uocomp.*sin(ocoarse);
    vcomp = urcomp.*sin(ocoarse) + uocomp.*cos(ocoarse);
% plot
    figure(2)
    hold on 
    quiver(xcoarse,ycoarse,ucomp,vcomp,'color',[0.4940 0.1840 0.5560],'LineWidth',1); %violet
    legend('Velocity Field');
    axis equal;
    xlabel('x');
    ylabel('y');
    xlim(dom);
    ylim(dom);
    %title('Velocity Field');
    hold off
%-- VELOCITY COMPONENTS -----------------------------------------------


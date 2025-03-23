clear all, close all

%-- PARAMETERS -------------------------------------------------------
% general
    mu=2.0;
% for domain
    dim = 101;
    dom = [-10,10];
    lessdim = 17;
%-- PARAMETERS -------------------------------------------------------


%-- SETTING DOMAIN ---------------------------------------------------
% general cartesian grid
    xgrid = linspace(dom(1),dom(2),dim)';
    ygrid = linspace(dom(1),dom(2),dim)';
    [xgrid,ygrid] = meshgrid(xgrid,ygrid);
% coarser cartesian grid for velocity
    xcoarse = linspace(dom(1),dom(2),lessdim)';
    ycoarse = linspace(dom(1),dom(2),lessdim)';
    [xcoarse,ycoarse] = meshgrid(xcoarse,ycoarse);
%-- SETTING DOMAIN ---------------------------------------------------


%-- STREAMLINES AND EQUIPOTENTIAL LINES ------------------------------
% functions
    phi = zeros(dim,dim);
    psi = zeros(dim,dim);
    for i=1:dim
        for j=1:dim
            if ( xgrid(i,j)^2 + ygrid(i,j)^2 ~= 0)
                phi(i,j) = + mu*xgrid(i,j) / (xgrid(i,j)^2 + ygrid(i,j)^2);
                psi(i,j) = - mu*ygrid(i,j) / (xgrid(i,j)^2 + ygrid(i,j)^2);
            else
                phi(i,j) = 0;
                psi(i,j) = 0;
            end
        end
    end
% const values for contour
    maxphi = 2;     
    maxpsi = 2;
    minphi = -2; 
    minpsi = -2; 
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
%-- STREAMLINES AND EQUIPOTENTIAL LINES ------------------------------


%-- VELOCITY COMPONENTS ----------------------------------------------
% functions
    ucomp = zeros(lessdim,lessdim);
    vcomp = zeros(lessdim,lessdim);
    for i=1:lessdim
        for j=1:lessdim
            if ( ycoarse(i,j)^2+xcoarse(i,j)^2 ~= 0 )
                ucomp(i,j) = + mu * ( ycoarse(i,j)^2-xcoarse(i,j)^2 )/(ycoarse(i,j)^2+xcoarse(i,j)^2)^2;
                vcomp(i,j) = - mu * (  2*ycoarse(i,j)*xcoarse(i,j)  )/(ycoarse(i,j)^2+xcoarse(i,j)^2)^2;
            else
                ucomp(i,j) = 0;
                vcomp(i,j) = 0;
            end
        end
    end
% plot
    figure(2)
    hold on 
    quiver(xcoarse,ycoarse,ucomp,vcomp,1,'color',[0.4940 0.1840 0.5560],'LineWidth',1); %violet
    legend('Velocity Field');
    axis equal;
    xlabel('x');
    ylabel('y');
    xlim(dom);
    ylim(dom);
    %title('Velocity Field');
    hold off
%-- VELOCITY COMPONENTS ----------------------------------------------


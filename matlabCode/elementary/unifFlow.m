clear all, close all

%-- PARAMETERS ---------------------------------------------------------
% general
    U = 20.0;  
% for domain
    dim = 101;
    dom = [-10,10];
    lessdim = 15;
%-- PARAMETERS ---------------------------------------------------------


%-- CARTESIAN DOMAIN ---------------------------------------------------
% general grid
    xgrid = linspace(dom(1),dom(2),dim)';
    ygrid = linspace(dom(1),dom(2),dim)';
    [xgrid,ygrid] = meshgrid(xgrid,ygrid);
% coarser grid for velocity
    xcoarse = linspace(dom(1),dom(2),lessdim)';
    ycoarse = linspace(dom(1),dom(2),lessdim)';
    [xcoarse,ycoarse] = meshgrid(xcoarse,ycoarse);
%-- CARTESIAN DOMAIN ---------------------------------------------------



%-- STREAMLINES AND EQUIPOTENTIAL LINES ---------------------------------
% functions
    phi = U*xgrid;
    psi = U*ygrid;
% const values for contour
    max1 = max(max(phi));
    max2 = max(max(psi));
    min1 = min(min(phi));
    min2 = min(min(psi));
    psivec = linspace(min1,max1,15)';
    phivec = linspace(min2,max2,15)';
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
%-- STREAMLINES AND EQUIPOTENTIAL LINES ---------------------------------


%-- VELOCITY COMPONENTS -------------------------------------------------
% functions
    ucomp = U*ones(lessdim,lessdim);
    vcomp = zeros(lessdim,lessdim);
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
%-- VELOCITY COMPONENTS -------------------------------------------------


clear all, close all

%-- PARAMETERS --------------------------------------------------------
% general
    b = 2.0;
% for domain
    dimr = 20;
    dimo = 40;
    domr = [b,b*8];
    dompPlot = [-domr(2)/2,domr(2)/2];
%-- PARAMETERS --------------------------------------------------------


%-- SETTING GRID ------------------------------------------------------
% general grid
    % set polar grid in CYLINDER space
    rgrid = linspace(domr(1),domr(2),dimr)';
    ogrid = linspace(0,2*pi,dimo)';
    [rgrid,ogrid] = meshgrid(rgrid,ogrid);
    % set cartesian grid in CYLINDER space
    egrid = rgrid.*cos(ogrid);
    ngrid = rgrid.*sin(ogrid);
% set cartesian grid in AIRFOIL space
    xxgrid = egrid + b^2 * egrid ./ (egrid.^2 + ngrid.^2);
    yygrid = ngrid - b^2 * ngrid ./ (egrid.^2 + ngrid.^2);
%-- SETTING GRID ------------------------------------------------------


%-- PLOT GRIDS --------------------------------------------------------
% uniform grid in e,n
    figure(1)
    hold on 
    axis equal;
    xlabel('\xi');
    ylabel('\eta');
    xlim(dompPlot);
    ylim(dompPlot);
    for i=1:size(egrid,1)
        plot(egrid(i,:),ngrid(i,:),'color',[0.4940 0.1840 0.5560],'LineWidth',1); %violet, horizontal 
    end
    for i=1:size(egrid,2) 
        plot(egrid(:,i),ngrid(:,i),'color',[0.4660 0.6740 0.1880],'LineWidth',1); %green, vertical 
    end
    hold off
% deformed grid in x,y
    figure(2)
    hold on 
    axis equal;
    xlabel('x');
    ylabel('y');
    xlim(dompPlot);
    ylim(dompPlot);
    for i=1:size(xxgrid,1)
        plot(xxgrid(i,:),yygrid(i,:),'color',[0.4940 0.1840 0.5560],'LineWidth',1); %violet, horizontal 
    end
    for i=1:size(xxgrid,2) 
        plot(xxgrid(:,i),yygrid(:,i),'color',[0.4660 0.6740 0.1880],'LineWidth',1); %green, vertical 
    end
    hold off
%-- PLOT GRIDS --------------------------------------------------------

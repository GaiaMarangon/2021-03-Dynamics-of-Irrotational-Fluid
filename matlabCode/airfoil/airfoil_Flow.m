clear all, close all


%-- PARAMETERS --------------------------------------------------------
% for shifting, rotating
    alpha = -pi/12;  % rotating
    beta = pi/18;   % y shifting
    e = 0.1;        % xshifting
    toll = 0.06;
% general
    U = -20.0;
    b = 2.0;
    a = b*(e+1);
    radius = a/cos(beta);
    gamma = 4*pi*radius*U*sin(alpha+beta);
    ncenter = a*tan(beta);
    ecenter = b*e;
% for grid dimensions
    dim = 301;
    ogni = 9;
    dom = [-6.5, 6.5];
%-- PARAMETERS --------------------------------------------------------


%-- SETTING GRID ------------------------------------------------------
% general grid
    % set cartesian grid in CYLINDER space
    egrid = linspace(-b*6,b*6,dim)';
    ngrid = linspace(-b*6,b*6,dim)';
    [egrid,ngrid] = meshgrid(egrid,ngrid);
    % set polar grid in CYLINDER space
    [rgrid,ogrid] = cart2pol(egrid,ngrid,ecenter,ncenter,alpha);
    % set to zero all points INTERNAL to cylinder
    [egrid,ngrid,rgrid,ogrid] = zeroInt(egrid,ngrid,rgrid,ogrid,radius,toll);
    % set cartesian grid in AIRFOIL space
    [xxgrid,yygrid] = cyl2air(egrid,ngrid,b);
    % set cartesian ROTATED grid in AIRFOIL space
    [xgrid,ygrid] = rotcart2cart(xxgrid,yygrid,alpha);
% special grid for plotting velocity
    % set polar adapted grid in CYLINDER space
    rPlot = linspace(radius,b*8,floor(dim/ogni));
    oPlot = linspace(0,2*pi,floor(dim/ogni));
    [rPlot,oPlot] = meshgrid(rPlot,oPlot);
    % set cartesian grid in CYLINDER space
    ePlot = rPlot.*cos(oPlot-alpha)+ecenter;       
    nPlot = rPlot.*sin(oPlot-alpha)+ncenter;        
    % set cartesian grid in AIRFOIL space
    [xxPlot,yyPlot] = cyl2air(ePlot,nPlot,b);
    % set cartesian ROTATED grid in AIRFOIL space
    [xPlot,yPlot] = rotcart2cart(xxPlot,yyPlot,alpha);
% data for plotting filled AIRFOIL
    % points of cylinder
    npoints = 100;
    interi=linspace(0,npoints,npoints)';
    theta=2*pi*interi/npoints;
    etheta = radius*cos(theta) + ecenter;
    ntheta = radius*sin(theta) + ncenter;
    % points of airfoil
    xxair = etheta + (b^2*etheta)./(etheta.^2 + ntheta.^2);
    yyair = ntheta - (b^2*ntheta)./(etheta.^2 + ntheta.^2);
    [xair,yair] = rotcart2cart(xxair,yyair,alpha);
%-- SETTING GRID ------------------------------------------------------
      
      
      
%-- CALL FOR FUNCTIONS  -----------------------------------------------
% for general space
    [phi,psi,phivec,psivec] = streamEqui(rgrid,ogrid,U,radius,gamma);
    [ucomp,vcomp,uair,vair] = velocity(rgrid,ogrid,egrid,ngrid,U,radius,gamma,b);
    [pcyl,pair,pveccyl,pvecair] = pressure(ucomp,vcomp,uair,vair,U);
% on surface
    [sucomp,svcomp,suair,svair] = velocity(radius*ones(npoints,1),theta+alpha,etheta,ntheta,U,radius,gamma,b);
    [spcyl,spair,spveccyl,spvecair] = pressure(sucomp,svcomp,suair,svair,U);
    [fairx,fairy] = forces(spair,xair,yair);
% reduced for plots
    [ucompPlot,vcompPlot,uairPlot,vairPlot] = velocity(rPlot,oPlot,ePlot,nPlot,U,radius,gamma,b);
%-- CALL FOR FUNCTIONS  -----------------------------------------------





%-- PLOTS FOR CYLINDER -----------------------------------------------
    figure(1)
    hold on 
    contour(egrid,ngrid,phi,phivec,'color',[0.0000 0.4470 0.7410],'LineWidth',1); %blue
    contour(egrid,ngrid,psi,psivec,'color',[0.4660 0.6740 0.1880],'LineWidth',1); %green
    fill(etheta,ntheta,'k','HandleVisibility','off')
    legend('Equipot. lines','Streamlines');
    axis equal;
    xlabel('\xi');
    ylabel('\eta');
    %title('Streamlines and Equipotential lines');
    hold off
% velocity
    figure(2)
    hold on 
    quiver(ePlot,nPlot,ucompPlot,vcompPlot,'color',[0.4940 0.1840 0.5560],'LineWidth',1); %violet
    fill(etheta,ntheta,'k','HandleVisibility','off')
    legend('Velocity Field');
    axis equal;
    xlabel('\xi');
    ylabel('\eta');
    %title('Velocity Field');
    hold off
% pressure
    figure(3)
    hold on 
    contourf(egrid,ngrid,pcyl,pveccyl,'LineColor','none');
    colormap(hsv)
    colorbar('westoutside');
    caxis([-2,max(pveccyl)])
    fill(etheta,ntheta,'k','HandleVisibility','off')
    legend('Pressure Field');
    axis equal;
    xlabel('\xi');
    ylabel('\eta');
    %title('Pressure Field');
    hold off
%-- PLOTS FOR CYLINDER -----------------------------------------------


%-- PLOTS FOR AIRFOIL -------------------------------------------------
% phi and psi
    figure(4)
    hold on 
    contour(xgrid,ygrid,phi,phivec,'color',[0.0000 0.4470 0.7410],'LineWidth',1); %blue
    contour(xgrid,ygrid,psi,psivec,'color',[0.4660 0.6740 0.1880],'LineWidth',1); %green
    fill(xair,yair,'k','HandleVisibility','off')
    legend('Equipot. lines','Streamlines');
    axis equal;
    xlabel('x');
    ylabel('y');
    xlim(dom)
    ylim(dom)
    %title('Streamlines and Equipotential lines');
    hold off
% velocity
    figure(5)
    hold on
    quiver(xPlot,yPlot,uairPlot,vairPlot,'color',[0.4940 0.1840 0.5560],'LineWidth',1); %violet
    fill(xair,yair,'k','HandleVisibility','off')
    legend('Velocity Field');
    axis equal;
    xlabel('x');
    ylabel('y');
    xlim(dom)
    ylim(dom)
    %title('Velocity Field');
    hold off
% pressure
    figure(6)
    hold on 
    contourf(xgrid,ygrid,pair,pvecair,'LineColor','none');
    colormap(hsv)
    colorbar('westoutside');
    caxis([-2,max(pvecair)])
    fill(xair,yair,'k','HandleVisibility','off')
    legend('Pressure Field');
    axis equal;
    xlabel('x');
    ylabel('y');
    xlim(dom)
    ylim(dom)
    %title('Pressure Field');
    hold off
% forces
    figure(7)
    hold on 
    fill(xair,yair,'k','HandleVisibility','off')
    quiver(xair,yair,fairx,fairy,'color',[0.3010 0.7450 0.9330],'LineWidth',1); %light blue
    legend('Force Field');
    axis equal;
    xlabel('x');
    ylabel('y');
    xlim(dom)
    ylim(dom)
    %title('Force Field');
    hold off
%-- PLOTS FOR AIRFOIL ------------------------------------------------
    
    
    
    
    


    
     
    
    
    
%-- FUNCTIONS --------------------------------------------------------
% from cartesian to adapted polar
function [rgrid,ogrid] = cart2pol(xgrid,ygrid,xcenter,ycenter,alpha)
    rgrid = sqrt( (xgrid-xcenter).^2 + (ygrid-ycenter).^2 );
    ogrid = atan2( (ygrid-ycenter),(xgrid-xcenter) ) + alpha;
end    
%----------------------------------------------------------------------    
% set to zero internal points
function [egrid,ngrid,rgrid,ogrid] = zeroInt(egrid,ngrid,rgrid,ogrid,radius,toll)
    dim = size(rgrid,1);
    for i=1:dim
        for j=1:dim
            if ( rgrid(i,j) < radius-toll)
                rgrid(i,j) = NaN;
                ogrid(i,j) = NaN;
                egrid(i,j) = NaN;
                ngrid(i,j) = NaN;
            end
        end
    end
end
%----------------------------------------------------------------------  
% from cylinder space to airfoil space with Joukowsky
function [xxgrid,yygrid] = cyl2air(egrid,ngrid,b)
    xxgrid = egrid + b^2 * egrid ./ (egrid.^2 + ngrid.^2);
    yygrid = ngrid - b^2 * ngrid ./ (egrid.^2 + ngrid.^2);
end   
%----------------------------------------------------------------------    
% from adapted (horizontal) cartesian to final cartesian
function [xgrid,ygrid] = rotcart2cart(xxgrid,yygrid,alpha)
    ss = ( xxgrid.^2 + yygrid.^2 ).^(1/2);
    tt = atan2(yygrid,xxgrid)+alpha;
    xgrid = ss.*cos(tt);
    ygrid = ss.*sin(tt);
end
%----------------------------------------------------------------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------
% streamfunction and equipotential function
function [phi,psi,phivec,psivec] = streamEqui(rgrid,ogrid,U,a,gamma)
    % function
    phi = U*( rgrid + a^2./rgrid ).*cos(ogrid) - gamma/(2*pi) .* ogrid;
    psi = U*( rgrid - a^2./rgrid ).*sin(ogrid) + gamma/(2*pi) .* log(rgrid/a);
    % const values for contour
    maxphi = 600;     
    maxpsi = 600;
    minphi = -600;
    minpsi = -600;
    phivec = linspace(minphi,maxphi,41)';
    psivec = linspace(minpsi,maxpsi,61)';
end
%---------------------------------------------------------------------------------------------------------
% velocity components
function [ucomp,vcomp,uair,vair] = velocity(r,o,e,n,U,a,gamma,b)
    %function in cylinder
    urcomp =  U*(1-a^2./r.^2).*cos(o);
    uocomp = -U*(1+a^2./r.^2).*sin(o) - gamma./(2*pi*r);
    ucomp  = urcomp.*cos(o) - uocomp.*sin(o);
    vcomp  = urcomp.*sin(o) + uocomp.*cos(o);
    % function in airfoil
    deno = (  e.^2 - n.^2 - b^2 ).^2 + 4*e.^2.*n.^2;
    A    = ( (e.^2 - n.^2).*(e.^2 - n.^2 - b^2) + 4*e.^2.*n.^2 ) ./ deno;
    B    =    - 2*e.*n*b^2 ./ deno;
    uair = ucomp.*A + vcomp.*B ;
    vair = vcomp.*A - ucomp.*B;
end
%---------------------------------------------------------------------------------------------------------
% (normalized) pressure function
function [pcyl,pair,pveccyl,pvecair] = pressure(ucomp,vcomp,uair,vair,U)
    % function
    pair = ( 1-(  uair.^2+vair.^2  )/U^2 )/2;
    pcyl = ( 1-( ucomp.^2+vcomp.^2 )/U^2 )/2;
    % const values for contour
    maxpcyl = +0.5;
    midpcyl = -0.9;
    minpcyl = -3.5;
    pveccyl = [linspace(minpcyl,midpcyl,150)';
              linspace(midpcyl,maxpcyl,15)']; 
    maxpair = +0.5;
    midpair = -0.9;
    minpair = -5;
    pvecair = [linspace(minpair,midpair,150)';
               linspace(midpair,maxpair,15)']; 
end
%---------------------------------------------------------------------------------------------------------
% forces
function [fairx,fairy] = forces(spair,xair,yair)
    % build shifted copies of vectors
    xair2 = xair;
    xair2(length(xair)) = [];
    xair(1) = [];
    yair2 = yair;
    yair2(length(yair)) = [];
    yair(1) = [];
    spair(length(spair)) = [];
    % use them to build force at each point
    ds = sqrt( (xair2-xair).^2 + (yair2-yair).^2 );
    orient = atan2( (xair2-xair), -(yair2-yair)  );
    fnorm = - spair .* ds;
    fairx = cos(orient).*fnorm;
    fairy = sin(orient).*fnorm;
    fairx = [fairx; 0];
    fairy = [fairy; 0];
end
%-- FUNCTIONS ---------------------------------------------------------

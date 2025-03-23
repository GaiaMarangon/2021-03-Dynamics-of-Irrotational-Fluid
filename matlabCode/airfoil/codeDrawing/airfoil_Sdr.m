clear all, close all


%-- PARAMETERS --------------------------------------------------
    b = 5.0;
    e = 0.1;        % x shifting
    beta = pi/18;   % y shifting
    a = b*(e+1);
    radius = a/cos(beta);
    ycenter = a*tan(beta);
    xcenter = b*e;
    alpha = pi/12;
%-- PARAMETERS --------------------------------------------------


%-- CYLINDER ----------------------------------------------------
% building points of cylinder
    npoints = 500;
    interi = linspace(0,npoints-1,npoints);
    theta = 2*pi*interi/(npoints-1);
    etheta = radius*cos(theta) + xcenter; 
    ntheta = radius*sin(theta) + ycenter; 
% building points for segments and angle
    p = [xcenter+radius*cos(pi/6); 
         ycenter+radius*sin(pi/6)];
    q = [xcenter+radius*cos(alpha); 
         ycenter-radius*sin(alpha)];
    ang = (pi/6+alpha)*linspace(0,20,21)/20;
    eang = 1.0*cos(ang-alpha) + xcenter; 
    nang = 1.0*sin(ang-alpha) + ycenter; 
% building points for flow
    U = 4;
    ux = -U*cos(alpha) * ones(4,1);
    uy =  U*sin(alpha) * ones(4,1);
    xu = ones(4,1)*12;
    yu = ycenter-tan(alpha)*(12-xcenter) * ones(4,1);
    yu(1) = yu(1)-3;
    yu(3) = yu(3)+3;
    yu(4) = yu(4)+6;
% plot cylinder
    figure(1)
    hold on
    axis equal
        % draw axis
        xdom=[-a-2,a+2+5];
        ydom=[-a-2,a+2];
        xaxx=linspace(xdom(1),xdom(2));
        yaxx=linspace(0,0);
        plot(xaxx,yaxx,'k-');
        xax=linspace(0,0);
        yax=linspace(ydom(1),ydom(2));
        plot(xax,yax,'k-');
        % circle    
        plot(etheta,ntheta,'r','LineWidth',2)
        plot(xcenter,ycenter,'xr','LineWidth',2)
        % segments and angle
        plot([xcenter,p(1)],[ycenter,p(2)],'-k','LineWidth',2);
        plot([xcenter,q(1)],[ycenter,q(2)],'--k','LineWidth',2);
        plot(eang,nang,'k','LineWidth',2);
        % flow
        quiver(xu,yu,ux,uy,0,'Color',[0 0.4470 0.7410],'LineWidth',2);
        % labels
        text(3.5,3.5,'r','FontWeight','bold','FontSize',16);
        text(2.1,1.3,'\theta','FontWeight','bold','FontSize',16);
        text(11,5,'U','FontWeight','bold','FontSize',16,'Color',[0 0.4470 0.7410]);
        text(11.7,-0.4,'\xi','FontSize',14)
        text(-0.8,7.2,'\eta','FontSize',14)
    xlim(xdom);
    ylim(ydom);
    hold off
%-- CYLINDER ----------------------------------------------------


%-- AIRFOIL 1 ---------------------------------------------------
% punti airfoil
    xair = etheta + (b^2*etheta)./(etheta.^2 + ntheta.^2);
    yair = ntheta - (b^2*ntheta)./(etheta.^2 + ntheta.^2);
% punti flow
    U = 6;
    ux = -U*cos(alpha) * ones(4,1);
    uy =  U*sin(alpha) * ones(4,1);
    xu = ones(4,1)*17.5;
    yu = [-8;-3;3;7];
% plot airfoil 
    figure(2)
    hold on
    axis equal
        % draw axis
        xdom=[-2.5*b,2.5*b+6];
        ydom=[-2.5*b,2.5*b-1];
        xaxx=linspace(xdom(1),xdom(2));
        yaxx=linspace(0,0);
        plot(xaxx,yaxx,'k-');
        xax=linspace(0,0);
        yax=linspace(xdom(1),xdom(2));
        plot(xax,yax,'k-');
        % airfoil
        plot(xair,yair,'r','LineWidth',2)
        % flow
        quiver(xu,yu,ux,uy,0,'Color',[0 0.4470 0.7410],'LineWidth',2);
        % labels
        text(15.5,8.5,'U','FontWeight','bold','FontSize',16,'Color',[0 0.4470 0.7410]);
        text(17,-0.9,'x','FontSize',14)
        text(-1.5,10.7,'y','FontSize',14)
    xlim(xdom);
    ylim(ydom);
    hold off
%-- AIRFOIL 1 ---------------------------------------------------

%-- AIRFOIL 2 ---------------------------------------------------
% punti airfoil
    [xxair,yyair] = rotcart2cart(xair,yair,alpha);
% punti flow
    ux = -ones(4,1)*6;
    uy = zeros(4,1);
    xu = ones(4,1)*17.5;
    yu = [-8;-3;2;7];
% punti assi secondari
    px = [ 17*cos(alpha), 17*sin(alpha)];
    qx = [-11*cos(alpha),-11*sin(alpha)];
    py = [-10*sin(alpha), 10*cos(alpha)];
    qy = [ 10*sin(alpha),-10*cos(alpha)];
% plot airfoil 
    figure(3)
    hold on
    axis equal
        % draw axis
        xdom=[-2.5*b,2.5*b+6];
        ydom=[-2.5*b,2.5*b-1];
        xaxx=linspace(xdom(1),xdom(2));
        yaxx=linspace(0,0);
        plot(xaxx,yaxx,'k-');
        xax=linspace(0,0);
        yax=linspace(xdom(1),xdom(2));
        plot(xax,yax,'k-');
        % airfoil
        plot(xxair,yyair,'r','LineWidth',2)
        % flow
        quiver(xu,yu,ux,uy,0,'Color',[0 0.4470 0.7410],'LineWidth',2);
        % secondary axis
        plot([px(1),qx(1)],[px(2),qx(2)],'--k');
        plot([py(1),qy(1)],[py(2),qy(2)],'--k');
        % labels
        text(15.5,8.5,'U','FontWeight','bold','FontSize',16,'Color',[0 0.4470 0.7410]);
        text(17,-0.9,'x''','FontSize',14)
        text(-1.5,10.7,'y''','FontSize',14)
        text(15,5,'x','FontSize',14)
        text(-3.7,9.2,'y','FontSize',14)
    xlim(xdom);
    ylim(ydom);
    hold off
%-- AIRFOIL 2 ---------------------------------------------------




%-- FUNCTIONS ---------------------------------------------------------    
% from adapted (rotated) cartesian to final cartesian
function [xgrid,ygrid] = rotcart2cart(xxgrid,yygrid,alpha)
    ss = ( xxgrid.^2 + yygrid.^2 ).^(1/2);
    tt = atan2(yygrid,xxgrid)+alpha;
    xgrid = ss.*cos(tt);
    ygrid = ss.*sin(tt);
end
%----------------------------------------------------------------------

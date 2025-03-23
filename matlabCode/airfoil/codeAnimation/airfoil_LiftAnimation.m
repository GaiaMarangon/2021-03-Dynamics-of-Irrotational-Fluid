clear all, close all

% FRAMES SETTING
% Define number of frames
nr_fr = 40;
% Initialize matrix using 'moviein'
frames = moviein(nr_fr); 
% Save in avi file
v1 = VideoWriter('airfoil_LiftAnimation.avi','Uncompressed AVI');
v1.FrameRate = 8;
open(v1)

%-- ITERATIONS FOR FRAMES - ANDATA -----------------------------------------------------------------
for k = 1 : nr_fr
    % CODE TO BE ITERATED
        %-- PARAMETERS --------------------------------------------------------
        % for shifting, rotating
            beta = pi/18;   % y shifting
            alpha = -beta + (pi/12+beta)*(k-1)/(nr_fr-1);
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
            dom = [-6.5, 6.5];
        %-- PARAMETERS --------------------------------------------------------


        %-- SETTING GRID ------------------------------------------------------
        % data for plotting filled AIRFOIL
            % points of cylinder
            npoints = 1000;
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
        % on surface
            [sucomp,svcomp,suair,svair] = velocity(radius*ones(npoints,1),theta+alpha,etheta,ntheta,U,radius,gamma,b);
            [spcyl,spair,spveccyl,spvecair] = pressure(sucomp,svcomp,suair,svair,U);
            [fairx,fairy] = forces(spair,xair,yair);
        %-- CALL FOR FUNCTIONS  -----------------------------------------------


        %-- LIFT -------------------------------------------------------------
        figure(1)
        clf
        hold on
        xlim([-1.5,1.5]);
        ylim([-1,12]);
            p = [0,0];
            for i=1:npoints
                quiver(p(1),p(2),fairx(i),fairy(i),0,'color',[0.3010 0.7450 0.9330],'LineWidth',1); %light blue
                p(1) = p(1)+fairx(i);
                p(2) = p(2)+fairy(i);
            end    
        quiver(0,0,p(1),p(2),0,'color',[0.0000 0.4470 0.7410],'LineWidth',1); 
        hold off
        %-- LIFT -------------------------------------------------------------
    % CODE TO BE ITERATED
    frame = getframe(gcf);
    writeVideo(v1,frame);
    if (k == 1 )     %first frame: keep same figure for some frames
        for h=1:10
            frame = getframe(gcf);
            writeVideo(v1,frame);
        end
    end
end
%-- ITERATIONS FOR FRAMES - ANDATA -----------------------------------------------------------------


%-- ITERATIONS FOR FRAMES - RITORNO -----------------------------------------------------------------
for k = 1 : nr_fr
    % CODE TO BE ITERATED 
        %-- PARAMETERS --------------------------------------------------------
        % for shifting, rotating
            beta = pi/18;   % y shifting
            alpha = -beta + (pi/12+beta)*(nr_fr-k)/(nr_fr-1);
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
            dom = [-6.5, 6.5];
        %-- PARAMETERS --------------------------------------------------------


        %-- SETTING GRID ------------------------------------------------------
        % data for plotting filled AIRFOIL
            % points of cylinder
            npoints = 1000;
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
        % on surface
            [sucomp,svcomp,suair,svair] = velocity(radius*ones(npoints,1),theta+alpha,etheta,ntheta,U,radius,gamma,b);
            [spcyl,spair,spveccyl,spvecair] = pressure(sucomp,svcomp,suair,svair,U);
            [fairx,fairy] = forces(spair,xair,yair);
        %-- CALL FOR FUNCTIONS  -----------------------------------------------


        %-- LIFT -------------------------------------------------------------
        figure(1)
        clf
        hold on
        xlim([-1.5,1.5]);
        ylim([-1,12]);
            p = [0,0];
            for i=1:npoints
                quiver(p(1),p(2),fairx(i),fairy(i),0,'color',[0.3010 0.7450 0.9330],'LineWidth',1); %light blue
                p(1) = p(1)+fairx(i);
                p(2) = p(2)+fairy(i);
            end    
        quiver(0,0,p(1),p(2),0,'color',[0.0000 0.4470 0.7410],'LineWidth',1); 
        hold off
        %-- LIFT -------------------------------------------------------------
    % CODE TO BE ITERATED
    frame = getframe(gcf);
    writeVideo(v1,frame);
    if (k == nr_fr )     %last frame: keep same figure for some frames
        for h=1:10
            frame = getframe(gcf);
            writeVideo(v1,frame);
        end
    end
end
%-- ITERATIONS FOR FRAMES - RITORNO -----------------------------------------------------------------
close(v1)






%-- FUNCTIONS --------------------------------------------------------    
% from adapted (horizontal) cartesian to final cartesian
function [xgrid,ygrid] = rotcart2cart(xxgrid,yygrid,alpha)
    ss = ( xxgrid.^2 + yygrid.^2 ).^(1/2);
    tt = atan2(yygrid,xxgrid)+alpha;
    xgrid = ss.*cos(tt);
    ygrid = ss.*sin(tt);
end
%----------------------------------------------------------------------
%----------------------------------------------------------------------
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

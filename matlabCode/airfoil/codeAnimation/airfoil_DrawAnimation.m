clear all, close all

% FRAMES SETTING
% Define number of frames
nr_fr = 20;
% Initialize matrix using 'moviein'
frames = moviein(nr_fr); 
% Save in avi file
v1 = VideoWriter('cylinder_ebetaAnimation.avi','Uncompressed AVI');
v2 = VideoWriter('airfoil_ebetaAnimation.avi','Uncompressed AVI');
v1.FrameRate = 8;
v2.FrameRate = 8;
open(v1)
open(v2)


%-- ITERATIONS FOR FRAMES - ANDATA -----------------------------------------------------------------
for k = 1 : nr_fr
    
    % CODE TO BE ITERATED
        % parametri
        b = 5.0;
        %e = 0;                    
        e = 0.1*(k-1)/(nr_fr-1);    %this one is varying
        a = b*(e+1);
        %beta = 0;               
        beta = pi/18*(k-1)/(nr_fr-1);   %this one is varying
        % punti del cerchio
        radius = a/cos(beta);
        ycenter = +a*tan(beta);
        xcenter = b*e;
        npoints = 500;
        interi = linspace(0,npoints-1,npoints);
        theta = 2*pi*interi/(npoints-1);
        etheta = radius*cos(theta) + xcenter; 
        ntheta = radius*sin(theta) + ycenter; 
        % punti airfoil
        xair = etheta + (b^2*etheta)./(etheta.^2 + ntheta.^2);
        yair = ntheta - (b^2*ntheta)./(etheta.^2 + ntheta.^2);
        % plot cylinder
        figure(1)
        clf
        hold on
        axis equal
        % draw axis
            xdom=[-a-2,a+2];
            ydom=[-a-2,a+2];
            xaxx=linspace(xdom(1),xdom(2));
            yaxx=linspace(0,0);
            plot(xaxx,yaxx,'k-');
            xax=linspace(0,0);
            yax=linspace(ydom(1),ydom(2));
            plot(xax,yax,'k-');
        plot(etheta,ntheta,'b','LineWidth',1)
        plot(xcenter,ycenter,'xb','LineWidth',1)
        xlabel('\xi');
        ylabel('\eta');
        xlim(xdom);
        ylim(ydom);
        hold off
        %------------------------------------
        frame = getframe(gcf);
        writeVideo(v1,frame);
        if (k == 1 )     %first frame: keep same figure for some frames
            for h=1:10
                frame = getframe(gcf);
                writeVideo(v1,frame);
            end
        end
        % plot airfoil 
        figure(2)
        clf
        hold on
        axis equal
        xlabel('x');
        ylabel('y');
        xlim([-2.5*b,2.5*b]);
        ylim([-2.5*b,2.5*b]);
        % draw axis
            xaxx=linspace(-2.5*b,2.5*b);
            yaxx=linspace(0,0);
            plot(xaxx,yaxx,'k-');
            xax=linspace(0,0);
            yax=linspace(-2.5*b,2.5*b);
            plot(xax,yax,'k-');
        %title('Airfoil as \beta vary');
        xlabel('x');
        ylabel('y');
        plot(xair,yair,'r','LineWidth',1);
        hold off
        %------------------------------------
        frame = getframe(gcf);
        writeVideo(v2,frame);
        if (k == 1 )     %first frame: keep same figure for some frames
            for h=1:10
                frame = getframe(gcf);
                writeVideo(v2,frame);
            end
        end
    % END CODE TO BE ITERATED
end
%-- ITERATIONS FOR FRAMES - ANDATA -----------------------------------------------------------------





%-- ITERATIONS FOR FRAMES - RITORNO -----------------------------------------------------------------
for k = 1 : nr_fr
    
   % CODE TO BE ITERATED
        % parametri
        b = 5.0;
        %e = 0;                    
        e = 0.1*(nr_fr-k)/(nr_fr-1);    %this one is varying
        a = b*(e+1);
        %beta = 0;               
        beta = pi/18*(nr_fr -k)/(nr_fr-1);   %this one is varying
        % punti del cerchio
        radius = a/cos(beta);
        ycenter = +a*tan(beta);
        xcenter = b*e;
        npoints = 500;
        interi = linspace(0,npoints-1,npoints);
        theta = 2*pi*interi/(npoints-1);
        etheta = radius*cos(theta) + xcenter; 
        ntheta = radius*sin(theta) + ycenter; 
        % punti airfoil
        xair = etheta + (b^2*etheta)./(etheta.^2 + ntheta.^2);
        yair = ntheta - (b^2*ntheta)./(etheta.^2 + ntheta.^2);
        % plot cylinder
        figure(1)
        clf
        hold on
        % draw axis
            xdom=[-a-2,a+2];
            ydom=[-a-2,a+2];
            xaxx=linspace(xdom(1),xdom(2));
            yaxx=linspace(0,0);
            plot(xaxx,yaxx,'k-');
            xax=linspace(0,0);
            yax=linspace(ydom(1),ydom(2));
            plot(xax,yax,'k-');
        plot(etheta,ntheta,'b','LineWidth',1)
        plot(xcenter,ycenter,'xb','LineWidth',1)
        axis equal
        xlabel('\xi');
        ylabel('\eta');
        xlim(xdom);
        ylim(ydom);
        hold off
        %------------------------------------
        frame = getframe(gcf);
        writeVideo(v1,frame);
        if (k == nr_fr )     %last frame: keep same figure for some frames
            for h=1:10
                frame = getframe(gcf);
                writeVideo(v1,frame);
            end
        end
        % plot airfoil 
        figure(2)
        clf
        hold on
        axis equal
        xlabel('x');
        ylabel('y');
        xlim([-2.5*b,2.5*b]);
        ylim([-2.5*b,2.5*b]);
        % draw axis
            xaxx=linspace(-2.5*b,2.5*b);
            yaxx=linspace(0,0);
            plot(xaxx,yaxx,'k-');
            xax=linspace(0,0);
            yax=linspace(-2.5*b,2.5*b);
            plot(xax,yax,'k-');
        %title('Airfoil as \beta vary');
        plot(xair,yair,'r','LineWidth',1);
        hold off
        %------------------------------------
        frame = getframe(gcf);
        writeVideo(v2,frame);
        if (k == nr_fr )     %last frame: keep same figure for some frames
            for h=1:10
                frame = getframe(gcf);
                writeVideo(v2,frame);
            end
        end
    % END CODE TO BE ITERATED
end
%-- ITERATIONS FOR FRAMES - RITORNO -----------------------------------------------------------------
close(v1)
close(v2)

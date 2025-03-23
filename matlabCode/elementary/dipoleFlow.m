clear all, close all

%-- PARAMETERS ----------------------------------------------------
%general
    a=2;
    m=20;
% for domain
    dim = 100;
    dom = [-10,10];
    lessdim = 21;
%-- PARAMETERS ----------------------------------------------------


%-- SETTING DOMAIN -------------------------------------------------
% general cartesian grid
    xgrid = linspace(dom(1),dom(2),dim)';
    ygrid = linspace(dom(1),dom(2),dim)';
    [xgrid,ygrid] = meshgrid(xgrid,ygrid);
% coarser cartesian grid for velocity
    xcoarse = linspace(dom(1),dom(2),lessdim)';
    ycoarse = linspace(dom(1),dom(2),lessdim)';
    [xcoarse,ycoarse] = meshgrid(xcoarse,ycoarse);
%-- SETTING DOMAIN -------------------------------------------------



%-- STREAMLINES AND EQUIPOTENTIAL LINES ----------------------------
% functions
    phi = zeros(dim,dim);
    psi = zeros(dim,dim);
    for i=1:dim
        for j=1:dim
            %phi
            if ( (xgrid(i,j)+a)^2 + ygrid(i,j)^2 >  0  && (xgrid(i,j)-a)^2 + ygrid(i,j)^2  > 0 )
                 phi(i,j) = m/(2*pi) * log( sqrt( (xgrid(i,j)+a)^2+ygrid(i,j)^2 )/sqrt( (xgrid(i,j)-a)^2+ygrid(i,j)^2 ) );
            elseif ( (xgrid(i,j)+a)^2 + ygrid(i,j)^2 ==  0  )
                phi(i,j) = -m*1.0e14;
            else
                phi(i,j) = m*1.0e14;
            end 
            %psi
            if ((xgrid(i,j) ~= -a) && (xgrid(i,j) ~= a))
                    psi(i,j) = m/(2*pi) * ( atan2( ygrid(i,j),(xgrid(i,j)+a) ) - atan2( ygrid(i,j),(xgrid(i,j)-a) ) );
            elseif (xgrid(i,j) == -a)
                    if (ygrid(i,j)>0)
                        psi(i,j) = m/(2*pi) * ( +pi/2 - atan2( ygrid(i,j),(xgrid(i,j)-a) )  );
                    elseif(ygrid(i)<0)
                        psi(i,j) = m/(2*pi) * ( -pi/2 - atan2( ygrid(i,j),(xgrid(i,j)-a) )  );
                    else
                        psi(i,j) = m/(2*pi) * (       - atan2( ygrid(i,j),(xgrid(i,j)-a) )  );
                    end
            elseif (xgrid(i,j) == a)
                    if (ygrid(i,j)>0)
                        psi(i,j) = m/(2*pi) * ( atan2( ygrid(i,j),(xgrid(i,j)+a) ) - pi/2 );
                    elseif (ygrid(i,j)<0)
                        psi(i,j) = m/(2*pi) * ( atan2( ygrid(i,j),(xgrid(i,j)+a) ) + pi/2 );
                    else
                        psi(i,j) = m/(2*pi) * ( atan2( ygrid(i,j),(xgrid(i,j)+a) )        );
                    end    
            end
        end
    end
% const values for contour
    maxphi = 5;  
    minphi = -5; 
    phivec = linspace(minphi,maxphi,25)';
% plot
    figure(1)
    hold on
    contour(xgrid,ygrid,phi,phivec,'color',[0.0000 0.4470 0.7410],'LineWidth',1); %blue
    contour(xgrid,ygrid,psi,35,'color',[0.4660 0.6740 0.1880],'LineWidth',1); %green
    legend('Equipot. lines','Streamlines');
    axis equal;
    xlabel('x');
    ylabel('y');
    xlim(dom);
    ylim(dom);
    %title('Streamlines and Equipotential lines');
    hold off
%-- STREAMLINES AND EQUIPOTENTIAL LINES ----------------------------



%-- VELOCITY COMPONENTS --------------------------------------------
% functions
    ucomp = m/(2*pi) * ( (xcoarse+a)./( (xcoarse+a).^2+ycoarse.^2 ) - (xcoarse-a)./( (xcoarse-a).^2+ycoarse.^2 ) );
    vcomp = m/(2*pi) * (    ycoarse ./( (xcoarse+a).^2+ycoarse.^2 ) -    ycoarse ./( (xcoarse-a).^2+ycoarse.^2 ) );
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
%-- VELOCITY COMPONENTS --------------------------------------------


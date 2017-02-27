function h = create_filled_ellipse(xc,yc,rx,ry,c,N)
% This function creates a filled ellipse. 
% xc,yc - center
% rx, ry - axis lengths along x and y
% c - color array
% N - number of points along boundary (default value is 100)

    if nargin < 6
            N = 100;
    end
    % create a set of angular coordinates
    theta = linspace(0,2*pi,N);
    
    % calculate x and y co-ordiantes on boundary
    x = xc + rx*cos(theta);
    y = yc + ry*sin(theta);

    % create a patch object with the specified color
    h = patch(x,y,c);
    
end

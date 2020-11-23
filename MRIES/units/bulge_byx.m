function fig = bulge_byx(x1, y1, x2, y2, b,direction_mark,parentaxes)
% draw a circle arc from P1 = (x1, y1) to P2 = (x2, y2) that bulges a 
% distance b to the LEFT of the arrow from P1 to P2. (use b < 0 to get a 
% bulge to the right.) The radius of the circle will be called "r". 

% Setting up how test appears for Matlab
set(0, 'DefaultTextInterpreter', 'latex', 'DefaultTextFontSize', 10)

if (b == 0)
    error('code does not handle the straight-line case');
end
P1 = [x1, y1]';           % The ' means "transpose", so P1 is a col. vector
P2 = [x2, y2]';
v = P2 - P1;              % vector from p1 to p2
t2 = sqrt(dot(v, v));     % distance between p1 and p2
if (t2 == 0)
    error('start and end point must be distinct')
end
v = v / t2;               % make v have unit length
vp = [-v(2) v(1)]';        % v rotated 90 degrees clockwise. 
M = (P1 + P2)/2;          % midpoint of chord. 
t = t2/2;                 % distance from M to P1 or to P2. 

% Let S denote the midpoint of the arc we're drawing. 
%
% The distance u from M to the circle center C satisfies two
% equations: 
% First, 
%   u^2 + t^2 = r^2  
% because C M P2 is a right triangle, and the distance from C to P2 is r,
% the radius of the circle
%
% Second
%   u + |b| = r
% because the union of the segments CM and MS is a radius of the circle,
% and CM has length u, while MS has length b. 
%
% We can solve for u:
% u = (t^2 - b^2)/(2|b|). 
%
u = (t^2 - b^2)/(2*abs(b)); 

% Now compute the locations of S and C, and the radius. 
S1 = M + b*vp;
S2 = M - b*vp;
if (sqrt(dot(S1,S1))>sqrt(dot(S2,S2)))
    b = -b;
end
S = M + b*vp;
if (b > 0)
    C = M - u*vp;
else
    C = M + u * vp;
end

r = sqrt(u^2 + t^2); 

% From here on, it's downhill sledding. 
% Find the ray from C to P1 and C to P2:

dir1 = (P1 - C)/r;
dir2 = (P2 - C)/r;

% Find the angles of these rays relative to the x-axis
% (so a ray in the positive-x direciton ahs angle 0; postive y gives pi/2
% negative y gives -pi/2, and so on, with the max values being +pi and -pi.
% 

theta1 = atan2(dir1(2), dir1(1));
theta2 = atan2(dir2(2), dir2(1));

% We now want to build points of the form 
% C + r [cos t, sin t]'
%
% where t ranges from theta1 to theta2. 
%
% But the wraparound at -pi and pi poses problems. 
%
% One observation: for b > 0, the angle t should be DECREASING from 
% theta1 to theta2, so if theta1 is less than theta2, we need to add 2pi. 
% The opposite holds if b < 0. 

if (b > 0) 
    if (theta1 < theta2)
        theta1 = theta1 + 2*pi;
    end
else
    if (theta2 < theta1)
        theta2 = theta2 + 2*pi; 
    end
end

n = 100; % draw 15 points. 
dt = (theta2 - theta1) / (n-1); 
colorvec = jet(n);
% if direction_mark == 1 means that it is bidirectional.
if direction_mark == 1
%     colorvec = zeros(n,3);   %bidirectional connection is black.
    colorvec = ones(n,1)*[0 0 0];      % bidirectional connection is red.
end
xs = zeros(1, n);
ys = zeros(1, n);
for i = 1:n  % "for i = 0 to n" in many programming languages
    u = theta1 + (i-1)*dt; 
    vv = (r * [cos(u), sin(u)]') + C;
    xs(i) = vv(1);
    ys(i) = vv(2);
end
for j = 1:n-1
      fig(j) = plot([xs(j) xs(j+1)],[ys(j) ys(j+1)],'color',colorvec(j,:),'linewidth',2,'Parent',parentaxes);

%     fig(j) = plot([xs(j) xs(j+1)],[ys(j) ys(j+1)],'color',colorvec(j,:),'linewidth',2);
end













function fig = plot_3Dcurve_colordir(p0,p1,linewidth,colorvec)
if nargin < 3
    linewidth = 2;
end

t = 0:0.1:1;
n = 1;
H=10;
for j = t
    a(n)=p0(1)+(p1(1)-p0(1))*j;
    b(n) = p0(2)+(p1(2)-p0(2))*j;
    h(n) = p0(3)+(p1(3)-p0(3))*j+H*sin(pi*j);
    n= n + 1;
end

if nargin < 4
    colorvec = colormap(jet(length(a)));
elseif nargin == 4
    colorvec = repmat(colorvec,[length(a),1]);
end

if length(linewidth) == 1
    linewidth = repmat(linewidth,[length(a),1]);
elseif length(linewidth) == 2
    lw = linewidth;
    linewidth = linspace(lw(1),lw(2),length(a));
end

for j = 1:length(a)-1
    
    fig(j,:) = line([a(j) a(j+1)],[b(j) b(j+1)],[h(j) h(j+1)],'Color',colorvec(j,:),'LineWidth',linewidth(j));
end

end
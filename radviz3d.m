function [U,D,theta,L] = radviz3d(D,mnv,mxv,Z,suppress_range_warnings)

% [U,D,theta,L] = radviz3d(D,mnv,mxv,Z,suppress_warnings)
%
% Generates 3d-radviz plot of data
%
% See:
% A. Ibrahim, S. Rahnamayan, M. V. Martin, and K. Deb, 
% 3D-RadVis: Visualization of Pareto front in many-objective optimization, 
% in Evolutionary Computation (CEC), 2016 IEEE Congress on, 2016, 
% pp. 736-745.
%
% D = data set to visualise, m columns of features (objectives)
% and n rows of patterns
%
% optional inputs
%
% mnv = minimum values possible for features (1 by m vector) -- useful to
%       pass in if you know what the acheivable minima are, and they are 
%       not present in D. If not passed, calculated from D.
% mxv = maximim values possible for features (1 by m vector)  -- useful to
%       pass in if you know what the acheivable maxima are, and they are 
%       not present in D. If not passed, calculated from D. 
% Z = OPTIONAL input, denotes guides (e.g. in a decomposition algorithm, to 
%       also be denoted in the visualisation). If not passed, empty set 
%       used 
% suppress_range_warnings = optional boolean flag to suppress warnings if
%       data in D exceeds the bounds passed in the optional mnv and mxv 
%       arguments. If not passed as an argument, set to true. 
%
%
% (c) Jonathan Fieldsend, University of Exeter, 2017, 2020


[n,m] = size(D);
% extract min and max from data if not passed in
if exist('mxv','var') == 0
    fprintf('Calculating maxima\n');
    mxv = max(D);
elseif sum((max(D)-eps)>mxv)>0 && exist('suppress_range_warnings','var')==0
    if suppress_range_warnings == false
        error('Array of maximum values contains at least one element smaller than held in D');
    end
end
if exist('mnv','var') == 0
    fprintf('Calulating minima\n');
    mnv = min(D);
elseif sum((min(D)+eps)<mnv)>0 && exist('suppress_range_warnings','var')==0
    if suppress_range_warnings == false
        error('Array of minimum values contains at least one element larger than held in D');
    end
end

%fix equality for divisions
I = find(mxv==mnv);
range = mxv-mnv;
range(range<eps) = 1;

% normalise data
D = (D-repmat(mnv,n,1))./repmat(range,n,1);
[U,theta] = map_to_radial_coordinates(D,n,m);

% get distance to simplex
L = sum(D,2)-1;

% plot radviz
scatter3(U(:,1),U(:,2),L,50,abs(L));
hold on
if exist('Z','var') % plot attractor lines
    if isempty(Z)==0
    size(Z)
    [n2,m] = size(Z);
    %Z = (Z-repmat(mnv,n2,1))./(repmat(mxv,n2,1)-repmat(mnv,n2,1));
    [Uz,thetaz] = map_to_radial_coordinates(Z,n2,m);
    %plot3(Uz(:,1),Uz(:,2),zeros(n2,1),'r.');
    hold on
    
%     %plot vertical lines through attractors
%     r = [min(min(L)) max(max(L)) ];
%     if r(1)>0
%         r(1)=0;
%     end
%     if r(2)<0
%         r(2)=0;
%     end
%     for i=1:n2
%         plot3([Uz(i,1) Uz(i,1)],[Uz(i,2) Uz(i,2)],r,'r:');
%     end
    plot3(Uz(:,1),Uz(:,2),zeros(n2,1),'*r','MarkerSize',2);
    end
end


%plot the simplex
plot(cos(theta),sin(theta),'kx'); 
plot(cos([theta theta(1)]),sin([theta theta(1)]),'k--'); 
%plot bounding simplices
if min(min(L))<0
    plot3(cos([theta theta(1)]),sin([theta theta(1)]),ones(1,m+1)*min(min(L)),'k:');
end
if max(max(L))>0
    plot3(cos([theta theta(1)]),sin([theta theta(1)]),ones(1,m+1)*max(max(L)),'k:');
end
%plot vertical points
r = [min(min(L)) max(max(L)) ];
if r(1)>0
    r(1)=0;
end
if r(2)<0
    r(2)=0;
end
for i=1:m
    plot3([cos(theta(i)) cos(theta(i))],[sin(theta(i)) sin(theta(i))],r,'k:');
end

%colour by distance from simplex
axis square;
colormap('copper');
colorbar

axis([-1 1 -1 1 min(L) max(L)]) 
%axis([-1 1 -1 1 -1 5])
end

function [U,theta] = map_to_radial_coordinates(D,n,m)

U = zeros(n,2);

theta = 0:1:m-1;
theta = theta*(2*pi/m);

theta_m = repmat(theta,n,1);
U(:,1) = sum(D.*cos(theta_m),2)./sum(D,2);
U(:,2) = sum(D.*sin(theta_m),2)./sum(D,2);


end

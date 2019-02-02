%{
This is the code which should be run to actually perform the basepoint
filtration. This code is currently formatted as a script, but if the user
 wishes it could be quickly converted to a function. First thing is to 
make sure the paths for javaplex and ripser are both in the current 
directory, or added. 
%}

import edu.stanford.math.plex4.*;
addpath('../../ripser');


%Initially, a point cloid is required for the plotting
surface = load('cat0.mat');
V = [surface.surface.X, surface.surface.Y, surface.surface.Z];
triv = surface.surface.TRIV;

%{
Next, make sure that there is a computationally feasible
number of points, on current version likely no more than 600-700. Here
 we use the built in maxmin selector from JavaPlex to do so.
%}
nn = 500;
maxmin_selector = api.Plex4.createMaxMinSelector(V,nn);
W = V(maxmin_selector.getLandmarkPoints()+1,:);


n = size(V);
n = n(1);

M = size(triv);
M = M(1);

% This next section is for computing the distance matrix of the metric
% space. For some spaces, this could be as simple as computing the 
% Euclidean distance matrix. A function is given at the bottom of the 
% code for doing so. Sometimes calculating the distance matrix can be
% more involved however. In this specific example, we wanted to calculate
% the geodesic distance along the surface of the cat. This is first done
% by finding all Euclidean distances between points and then calculating an
% approximation of geodesic distance using djikstra's method when 
% only considering a few of the nearest neighbors to each point as 
% possessing an edge between them.

dW = zeros(nn);
%dWg = zeros(nn);
rowmin = zeros(nn,1);
for ii=1:(nn-1)
    for jj=(ii+1):nn
         dW(ii,jj) = distance(ii, jj, W);    
    end
end

dW = dW + dW';

indices = zeros(nn,1);

for ii=1:nn
    [tf, index]=ismember(W(ii,:),V,'rows');
    indices(ii) = index;
end

for ii=1:nn
    y = sort(dW(ii,:));
    rowmin(ii) = y(7);
end

clear B;
clear C;
for jj=1:nn
    B(jj,:) = dW(jj,:) > rowmin(jj);
end
dW(B) = 0;

D = dijkstra( dW , 1:nn );
D = (D + D') .*.5;

%{ 
This is the code which plots the points and allows the user to click
on a point which is then used as the basepoint for the filtration.
Note currently this is only supported for 3-dimensional point clouds. The 
user must ensure this, either by adding a 0 coordinate to each point for 
a 2-dimensional point cloud or by performing a dimension-reduction to 
the original point cloud through something such as PCA.

W' is point cloud and D is distance matrix
%}
basepointclick(W', D);


%this computes the Euclidean distance in R3
function edist = distance(index1,index2,points)
    x = points(index1,:);
    y = points(index2,:);
    edist = sqrt((x(1)-y(1))^2 + (x(2)-y(2))^2 + (x(3)-y(3))^2);
end
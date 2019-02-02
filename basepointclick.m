function h = basepointclick(pointCloud, dx)

if size(pointCloud, 1)~=3
    error('Input point cloud must be a 3*N matrix.');
end

% show the point cloud
h = gcf;
scatter3(pointCloud(1,:), pointCloud(2,:), pointCloud(3,:),20,'filled'); 
cameratoolbar('Show'); % show the camera toolbar
hold on; % so we can highlight clicked points without clearing the figure

% set the callback, pass pointCloud to the callback function
set(h, 'WindowButtonDownFcn', {@ripsercallback, pointCloud, dx}); 

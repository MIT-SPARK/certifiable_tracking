%% Example visualize script

% Pick category and instance
cat = "aeroplane";
inst = 4;

% Load
c = load(cat + '.mat', cat);
c = getfield(c,cat);
faces = c(inst).faces;
vertices = c(inst).vertices;

% Visualize
figure
trimesh(faces, vertices(:,1),vertices(:,2),vertices(:,3),'FaceAlpha',0.5)
axis equal
hold on
% for p = c(inst).pnames
%     coord = getfield(c(inst),p{1});
%     plot3(coord(1),coord(2),coord(3),'.r','MarkerSize',20)
% end
colormap(winter)
axis off

% plot line to keypoint
p = c(inst).pnames(5);
coord = getfield(c(inst),p{1});
centroid = mean(vertices);
plot3(coord(1),coord(2),coord(3),'.r','MarkerSize',20)
% quiver3(centroid(1),centroid(2),centroid(3),coord(1)-centroid(1),coord(2)-centroid(2),coord(3)-centroid(3),0,'LineWidth',3)
plot3([centroid(1); coord(1)], [centroid(2); coord(2)], [centroid(3); coord(3)],'LineWidth',3)
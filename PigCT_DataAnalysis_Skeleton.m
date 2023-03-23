%% BMEN E3820: Biomedical Engineering Laboratory II (Pig CT Data Analysis)
%  Written by: Dr. Lauren N. Heckelman
%  Spring 2023

%% Initialize the Workspace:
clear; clc; close all;

%% Load & View DICOM Images:
    Slash = '/'; % should be '/' on Mac, '\' on Windows
    files = dir(strcat(pwd, Slash, '*.dcm'));

    filenames = cell(1,length(files));

    for n = 1:length(filenames)
        filenames{1,n} = files(n).name;
    end

    filenames = natsort(filenames);

    testIM = dicomread(strcat(pwd, Slash, filenames{1,round(length(filenames)./2)}));
    thresh = 1250;
    IM = zeros(size(testIM,1), size(testIM,2), length(filenames));
    IM_mask = zeros(size(testIM,1), size(testIM,2), length(filenames));
    IM_overlay = zeros(size(testIM,1), size(testIM,2), length(filenames), 3);
    boundary_image= zeros(size(testIM,1), size(testIM,2), length(filenames));
    patella_whole = zeros(size(testIM,1), size(testIM,2), length(filenames));
    patella_bounary = zeros(size(testIM,1), size(testIM,2), length(filenames));
%     IM_corrected = 

    boundaries = cell(1,321); %321 is the size of Z.
    for n = 1:length(filenames)
        IM(:,:,n) = dicomread(strcat(pwd, Slash, filenames{1,n}));
        IM_mask(:,:,n) = IM(:,:,n) > thresh;
        IM_overlay(:,:,n,:) = imoverlay(IM(:,:,n)./3000,IM_mask(:,:,n), 'green');
        B = bwboundaries(IM_mask(:,:,n), 8, 'noholes');
        patella = IM_mask(:,:,n);
        [labeledImage, ~] = bwlabel(patella, 4);
        prop = regionprops(labeledImage, 'BoundingBox');
        props = cat(1, prop.BoundingBox);
        [~,ind] = maxk(props(:,3),2);
        second_max_ind = ind(end);
        patella(labeledImage ~= second_max_ind) = 0;
        patella_whole(:,:,n) = patella;
        Bpatella = bwboundaries(patella, 8, 'noholes');
        %[L,Centers] = imsegkmeans(IM_mask,2);
        %clustered_image(:,:,n) = L;
        for k = 1:length(B)
          boundary = B{k};
          for i = 1:length(boundary)
              boundary_image(boundary(i,1), boundary(i,2), n) = 1;
          end
          boundaries{n} = [[boundary(:,1)'; boundary(:,2)'], boundaries{n}];
        end
        %boundaries{n} = [B(:,2); B(:,1)]; %x| y
%         boundary_n = [];
%         for k = 1:length(B)
%            boundary = B{k};
%            x = [x, boundary(:,2)'];
%            y = [y, boundary(:,2)'];
%            % plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
%         end
%         boundary_n = [x ; y];
%         boundariess(n) = mat2cell(boundary_n);
    end
 

    numslices = size(IM,3);
       
    figure(1); clf;
    sliceViewer(IM, 'SliceNumber', 1);

    figure(2); clf;
    sliceViewer(IM_mask, 'SliceNumber', 1);

    figure(3); clf;
    sliceViewer(IM_overlay, 'SliceNumber', 3);

    figure(4); clf;
    sliceViewer(boundary_image,'SliceNumber', 1 )

    figure(5);
    slice = 50:140;
    patella_whole = patella_whole(:,:,slice);
    sliceViewer(patella_whole, 'SliceNumber', 1)
%%
% Use the pixel dimensions listed in the DICOM metadata to convert the 
% pixel index locations into physical coordinates in x, y, and z.

p_x = 1:512;
p_y = p_x;
p_z = 1:321;

res_y = 0.375;
res_x = res_y;
res_z = 0.6; % 1.8691588785 mm/px

[mesh_x, mesh_y, mesh_z] = meshgrid(p_x*res_x, p_y*res_y, p_z*res_z);

X = mesh_x.*IM_mask;
Y = mesh_y.*IM_mask; 
Z  = mesh_z.*IM_mask; 

figure(5)
scatter3(X(:), Y(:), Z(:) ,'.', 'MarkerEdgeColor', [227/255, 218/255, 201/255])
axis equal

%%
p_x = 1:512;
p_y = p_x;
p_z = 1:91;

res_y = 0.375;
res_x = res_y;
res_z = 0.6; % 1.8691588785 mm/px

[mesh_x, mesh_y, mesh_z] = meshgrid(p_x*res_x, p_y*res_y, p_z*res_z);

x = mesh_x.*patella_whole;
y = mesh_y.*patella_whole; 
z  = mesh_z.*patella_whole; 
shp = alphaShape(x(:), y(:), z(:));
shp_points = shp.Points;
xp = shp_points(:,1);
yp = shp_points(:,2);
zp = shp_points(:,3);

points = [xp, yp, zp];
writematrix( points,"patella_points.txt", 'delimiter', '\t')
[T, P] =  alphaTriangulation(shp);
stlwrite(triangulation(T,P),['foo.stl'])
figure(22)
plot(shp)
axis equal

% Write the 3D model to an STL file

[xp,yp,zp] = shp.Points;

% Convert the 3D model to triangles
tri = delaunay(xp(:),yp(:),zp(:));

% Write the 3D model to an STL file
stlwrite('sphere.stl',tri);
%stlwrite('test.stl',tri);

%make_STL_of_Array('lab7_2.stl', patella_whole, res_x, res_y, res_z);
%% Bone to print
z_slice = 50:180;
y_slice = 70:300;
x_slice = 70:178;


X_slice = X(x_slice, y_slice, z_slice);
Y_slice = Y(x_slice, y_slice, z_slice);
Z_slice = Z(x_slice, y_slice, z_slice);
IM_slice = IM_mask(x_slice, y_slice, z_slice);

%make_STL_of_Array('lab7_2.stl', IM_mask(x_slice, :, z_slice), res_x, res_y, res_z);


%%
figure(6)

scatter3(X_slice(:), Y_slice(:), Z_slice(:) ,'.', 'MarkerEdgeColor', [227/255, 218/255, 201/255])

%%

% n = 100;
% boundary = bwboundaries(IM_mask(:,:,n));
% IM(:,:,n) = dicomread(strcat(pwd, Slash, filenames{1,n}));
% IM_mask(:,:,n) = IM(:,:,n) > thresh;
% IM_overlay(:,:,n,:) = imoverlay(IM(:,:,n),IM_mask(:,:,n), 'green');
% 
% figure(4); clf;
% imshow(IM(:,:,n))
% plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
% \

%%
testing_img = IM_mask(:,:, 70);
[labeledImage, numberOfRegions] = bwlabel(testing_img, 4);
prop = regionprops(labeledImage, 'BoundingBox');
props = cat(1, prop.BoundingBox);
[mx,ind] = maxk(props(:,3),2);
second_max_ind = ind(2);
len = length(prop);
testing_img(labeledImage ~= second_max_ind) = 0;
figure(22);clf;
imshow(testing_img)
%%
% prop = regionprops(labeledImage, 'BoundingBox');
% props = cat(1, prop.BoundingBox);
% [row, col, v] = find(props == min(props(:,2)));
% props(row, :) = [];
% 
% 
% figure(18); clf;
% nexttile()
% imshow(testing_img)
% hold on
% for i = props'
%     xstart = floor(i(1));
%     ystart = floor(i(2));
%     xstep  = floor(i(3));
%     ystep = floor(i(4));
%     rectangle("Position", i, "EdgeColor", 'r')
%     testing_img(xstart:xstart+xstep, ystart:ystart+ystep) = 0;
% end
% hold off
% % rectangle("Position", [206.500000000000,198.500000000000,115,122], "EdgeColor", 'r');
% % hold on
% % rectangle("Position", [284.500000000000,97.5000000000000,32,34], "EdgeColor", 'b');
% % 
% % testing_img(206:206+115, 198:198+122) = 0;
% 
% nexttile()
% imshow(testing_img)
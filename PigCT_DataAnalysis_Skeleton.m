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
        IM_overlay(:,:,n,:) = imoverlay(IM(:,:,n),IM_mask(:,:,n), 'green');
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
        for k = 1:length(B)
          boundary = B{k};
          for i = 1:length(boundary)
              boundary_image(boundary(i,1), boundary(i,2), n) = 1;
          end
          boundaries{n} = [[boundary(:,1)'; boundary(:,2)'], boundaries{n}];
        end

    end
 

    numslices = size(IM,3);
       
    figure(1); clf;
    sliceViewer(IM, 'SliceNumber', 1);

    figure(2); clf;
    sliceViewer(IM_mask, 'SliceNumber', 1);

    figure(3); clf;
    sliceViewer(IM_overlay, 'SliceNumber', 3);

    figure(4); clf;
    sliceViewer(boundary_image,'SliceNumber', 1 );

    figure(5);
    slice = 50:140;
    patella_whole = patella_whole(:,:,slice);
    sliceViewer(patella_whole, 'SliceNumber', 1);
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

figure(6)
scatter3(X(:), Y(:), Z(:) ,'.', 'MarkerEdgeColor', [227/255, 218/255, 201/255])
title("3d point Cloud of Bone")
axis equal

%%


xp = mesh_x(:,:,slice).*patella_whole;
yp = mesh_y(:,:,slice).*patella_whole; 
zp  = mesh_z(:,:,slice).*patella_whole; 

figure(7)
scatter3(xp(:), yp(:), zp(:) ,'.', 'MarkerEdgeColor', [227/255, 218/255, 201/255])
title("3d Point Cloud of Patella")
axis equal

shp = alphaShape(xp(:), yp(:), zp(:));

shp_points = shp.Points;
xshp = shp_points(:,1);
yshp = shp_points(:,2);
zshp = shp_points(:,3);

points = [xshp, yshp, zshp];
writematrix( points,"patella_points.txt", 'delimiter', '\t')



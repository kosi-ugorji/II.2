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
%     IM_corrected = 
    boundariess = {};
    x = [];
    y = [];
    for n = 1:length(filenames)
        IM(:,:,n) = dicomread(strcat(pwd, Slash, filenames{1,n}));
        IM_mask(:,:,n) = IM(:,:,n) > thresh;
        IM_overlay(:,:,n,:) = imoverlay(IM(:,:,n),IM_mask(:,:,n), 'green');
        B = bwboundaries(IM_mask(:,:,n));
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
    sliceViewer(IM_overlay, 'SliceNumber', 1);

%%
% Use the pixel dimensions listed in the DICOM metadata to convert the 
% pixel index locations into physical coordinates in x, y, and z.

p_x = 1:512;
p_y = p_x;
p_z = 1:321;

res_y = 0.375;
res_x = res_y;
res_z = 0.6/ 321; % 1.8691588785 mm/px

[mesh_x, mesh_y] = meshgrid(p_x*res_x, p_y*res_y);
mesh_plot = mesh .*IM_mask(1:end, 1:end, 160);
figure(10)
scatter3(mesh_x, mesh_y, IM(1:end, 1:end, 10).*IM_mask(1:end, 1:end, 10)*res_z ,'.')
%%

n = 100;
boundary = bwboundaries(IM_mask(:,:,n));
IM(:,:,n) = dicomread(strcat(pwd, Slash, filenames{1,n}));
IM_mask(:,:,n) = IM(:,:,n) > thresh;
IM_overlay(:,:,n,:) = imoverlay(IM(:,:,n),IM_mask(:,:,n), 'green');

figure(4); clf;
imshow(IM(:,:,n))
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)


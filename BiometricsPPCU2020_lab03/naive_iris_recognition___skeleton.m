% iris segmentation, iris code, matching with Hamming distance
% initial file
% 2020. 11. 20.
% on the basis of Libor Masek's thesis: 
% https://www.peterkovesi.com/studentprojects/libor/index.html
% images: selected from the CASIA db, Chinese Academy of Sciences
% http://www.cbsr.ia.ac.cn/english/IrisDatabase.asp

close all;
clear;
clc;

filenames = {'./input/S1007R05.jpg', './input/S1007R07.jpg', ...
    './input/S1008L02.jpg', './input/S1008L10.jpg', ...
    './input/S1011L02.jpg', './input/S1011L09.jpg'};

I_templates = cell(0);
I_noise_mask = cell(0);


for idx = 1:length(filenames)
    iris = im2double(imread(filenames{idx}));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 1: calculate edge-images to circular Hough + linear Hough:
    iris_edges_all = imgaussfilt(iris, 9);
    iris_edges_all = edge(iris_edges_all, 'Canny');
    iris_horizontal_edges = edge(iris, 'Prewitt');
    
    
    
    
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 2: circular and linear Hough transforms to detect boundaries:
    [centers_1, radii_1] = imfindcircles(iris_edges_all,[30 70],'Method','TwoStage','Sensitivity',0.95);
    c_1 = centers_1(1,:);
    r_1 = radii_1(1);
    
    [centers_2, radii_2] = imfindcircles(iris_edges_all,[70 125],'Method','PhaseCode','Sensitivity',0.98);
    c_2 = centers_2(1,:);
    r_2 = radii_2(1);
    
    lower_part = round(3*size(iris,1)/4); % this name is used in the printer part (exc4)
    bottom = iris_horizontal_edges(offset:end,:);
    [H,T,R] = hough(bottom,'Theta',[-90:-85,85:89]);
    p = houghpeaks(H,5);
    H_lines_bottom = houghlines(bottom,T,R,p);
    H_bottom_row_coords = [];
    for k = 1:length(H_lines_bottom)
        H_bottom_row_coords = [H_bottom_row_coords, H_lines_bottom(k).point1(2), H_lines_bottom(k).point2(2)];
    end
    highest_row_to_remove = sum(offset,min(H_bottom_row_coords));
    
    offset = round(size(iris,1)/4);
    top = iris(1:offset,:);
    [H,T,R] = hough(top,'Theta',[-90:85,85:89]);
    p = houghpeaks(H,5);
    H_lines_top = houghlines(top,T,R,p);
    H_top_row_coords = [];
    for k = 1:length(H_lines_top)
        H_top_row_coords = [H_top_row_coords, H_lines_top(k).point1(2), H_lines_top(k).point2(2)];
    end
    lowest_row_to_remove = min(H_top_row_coords);
    lower_part = min(H_top_row_coords);
    masked_iris = iris;
    masked_iris(1:lowest_row_to_remove,:) = NaN;
    masked_iris(highest_row_to_remove:end,:) = NaN;
    
    
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 3: normalisation step
    % on the basis of Daugman's rubber sheet model
    
    
    
    
    
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 4: visualize:
    figure(1);
    subplot(length(filenames), 6, (idx-1)*6+1);
    imshow(iris);
    
    
    % please uncomment if you are ready with STEP 1
    subplot(length(filenames), 6, (idx-1)*6+2);
    imshow(iris_edges_all); % result of canny
    subplot(length(filenames), 6, (idx-1)*6+3);
    imshow(iris_horizontal_edges);
    
    
    
    % please uncomment if you are ready with STEP 2 / circular Hough
    subplot(length(filenames), 6, (idx-1)*6+4);
    imshow(iris);
    hold on;
    viscircles(c_1(1, :), r_1(1));
    viscircles(c_2(1, :), r_2(1));
    
    
    % please uncomment if you are ready with STEP 2 / linear Hough
    for k = 1:length(H_lines_bottom)
        plot([H_lines_bottom(k).point1(1), H_lines_bottom(k).point2(1)], ...
            [H_lines_bottom(k).point1(2), H_lines_bottom(k).point2(2)]+lower_part, 'g');
    end
    for k = 1:length(H_lines_top)
        plot([H_lines_top(k).point1(1), H_lines_top(k).point2(1)], ...
            [H_lines_top(k).point1(2), H_lines_top(k).point2(2)], 'b');
    end
    
    
    %{
    % please uncomment if you are ready with STEP 3
    subplot(length(filenames), 6, (idx-1)*6+5);
    imshow(polar_array);
    subplot(length(filenames), 6, idx*6);
    imshow(noise_mask);
    %}
    
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 5: calculate an iris-code - log Gabor
    
    
    
    
    
    
    
    
    
    
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 6: matching with Hamming distance














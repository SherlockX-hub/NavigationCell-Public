%% markPartSFP_color_v3.m
% This code is used to plot contour of part cells' spatial footprints.
% Trim or release version of v2.
% Input:
%       SFP_all: m*n*p matrix, SFP of all cells;
%       SFP_mark: m*n*q matrix (q < p), SFP of colored cells;
%       color: 1*3 matrix, [R, G, B].

% Created by Xiang Zhang, 2021.

function f = markPartSFP_color_v3(SFP_all,SFP_mark,color)
    
    if isempty(color)
        R = 1;
        G = 0.1;
        B = 0.1;
    elseif length(color) == 3
        R = color(1);
        G = color(2);
        B = color(3);
    else
        error('wrong input with color');
    end
    color_background = [1 1 1];
    color_blankCell = [225,232,240] / 255;
    color_cellContour = [116,112,129] / 255;
    
    SFP_all_temp = max(permute(SFP_all,[2 3 1]),[],3);
    SFP_mark_temp = max(permute(SFP_mark,[2 3 1]),[],3);
    
    % filled the specific cells with one color;
    SFP_all_temp(SFP_all_temp ~= 0) = 1; % all cells;
    SFP_mark_temp(SFP_mark_temp ~= 0) = 1; % marked cells;
    
    SFP_all_new_first = zeros(size(SFP_all_temp));
    SFP_all_new_second = zeros(size(SFP_all_temp));
    SFP_all_new_third = zeros(size(SFP_all_temp));
    SFP_all_new_label = SFP_all_temp; % 1 means the blank cells; use new variable to avoid the same color number;
    
    SFP_all_new_first(SFP_mark_temp == 1) = R; % marked cells;
    SFP_all_new_second(SFP_mark_temp == 1) = G;
    SFP_all_new_third(SFP_mark_temp == 1) = B;
    SFP_all_new_label(SFP_mark_temp == 1) = 0.5; % 0.5 means the marked cells;
    
    SFP_all_new_first(SFP_all_new_label == 1) = color_blankCell(1); % blank cells;
    SFP_all_new_second(SFP_all_new_label == 1) = color_blankCell(2);
    SFP_all_new_third(SFP_all_new_label == 1) = color_blankCell(3);
    
    SFP_all_new_first(SFP_all_new_label == 0) = color_background(1); % 0 means background;
    SFP_all_new_second(SFP_all_new_label == 0) = color_background(2);
    SFP_all_new_third(SFP_all_new_label == 0) = color_background(3);
    
    SFP_all_new(:,:,1) = SFP_all_new_first;
    SFP_all_new(:,:,2) = SFP_all_new_second;
    SFP_all_new(:,:,3) = SFP_all_new_third;
    
    f = figure;
    imagesc(SFP_all_new);
    
    hold on;
    for i = 1:size(SFP_all,1)
        SFP_all_temp_i = SFP_all(i,:,:);
        SFP_all_temp_i = permute(SFP_all_temp_i,[2 3 1]);
        [~,n,o] = size(SFP_all);
        SFP_all_contour = zeros(n,o);
        %     SFP_mark_contour = zeros(size(SFP_all,2:3));
        SFP_all_contour(SFP_all_temp_i ~= 0) = 1;
        contour(SFP_all_contour,[0.01 0.01],'color',color_cellContour,'LineWidth',1);
    end
    hold off;
    
    axis equal;
    axis off;
    
end
%% markMultiSFP_color_v3.m
% This code is used to plot SFP of cells with multiple identities.

% Input:
%        SFP_all: a m*n*p matrix, the SFP of all cells;
%        SFP_mark: a l*1 cell, each cell contains the SFP of each cell type;
%        SFP_name: a l*1 cell, each cell contains the name of each cell type;
%        SFP_color: a l*1 cell, each cell contains the color of each cell type;
%        islegend: a bool, whether show legend or not;

% v3: 2021.
% This version draw multiple identities in each cell.
% There are some logic mistakes in this codes.

% Created by Xiang Zhang, 2021.

function markMultiSFP_color_v3(SFP_all,SFP_mark,SFP_name,SFP_color,islegend)
    % SFP_mark, SFP_name, color are "cell" type, SFP_all are "matrix" type.
    
    if ~iscell(SFP_mark) || ~iscell(SFP_color) % || ~isstring(SFP_name)
        error("Error input, please check it again.");
    end
    
    if length(SFP_mark) ~= length(SFP_color) || length(SFP_mark) ~= length(SFP_name)
        error("Error length of input, please check it again.");
    end
    
    color_background = [1 1 1];
    color_blankCell = [225,232,240] / 255;
    color_cellContour = [116,112,129] / 255;
    
    SFP_all_temp = max(permute(SFP_all,[2 3 1]),[],3);
    SFP_all_temp(SFP_all_temp ~= 0) = 1; % all cells;
    
    SFP_all_new_first = zeros(size(SFP_all_temp));
    SFP_all_new_second = zeros(size(SFP_all_temp));
    SFP_all_new_third = zeros(size(SFP_all_temp));
    SFP_all_new_label = SFP_all_temp; % 1 means the blank cells; use new variable to avoid the same color number;
    
    for mark_i = 1:length(SFP_mark)
        SFP_mark_temp = SFP_mark{mark_i};
        SFP_color_temp = SFP_color{mark_i};
        SFP_color_num = size(SFP_color_temp,1);
        
        % only one color;
        if SFP_color_num == 1
            SFP_mark_temp = max(permute(SFP_mark{mark_i},[2 3 1]),[],3);
            SFP_mark_temp(SFP_mark_temp ~= 0) = 1; % marked cells;
            
            SFP_all_new_first(SFP_mark_temp == 1) = SFP_color{mark_i}(1); % marked cells;
            SFP_all_new_second(SFP_mark_temp == 1) = SFP_color{mark_i}(2);
            SFP_all_new_third(SFP_mark_temp == 1) = SFP_color{mark_i}(3);
            SFP_all_new_label(SFP_mark_temp == 1) = 0.5; % 0.5 means the marked cells;
            
        else % more than one color;
            
            for cell_i = 1:size(SFP_mark_temp,1)
                SFP_cell_temp = permute(SFP_mark_temp(cell_i,:,:), [2 3 1]);
                [SFP_cell_row, SFP_cell_col] = find(SFP_cell_temp ~= 0); % marked cells;
                SFP_cell_col_min = min(SFP_cell_col);
                SFP_cell_col_max = max(SFP_cell_col);
                
                SFP_cell_col_inter = ceil((SFP_cell_col_max - SFP_cell_col_min) / SFP_color_num);
                if SFP_cell_col_inter == 0, SFP_cell_col_inter = 1; end
                for col_i = 1:SFP_color_num - 1
                    SFP_cell_col_idx = find(SFP_cell_col >= SFP_cell_col_min + SFP_cell_col_inter * (col_i - 1) ...
                        & SFP_cell_col < SFP_cell_col_min + SFP_cell_col_inter * (col_i));
                    
                    % marked cells;
                    for SFP_cell_col_idx_i = 1:length(SFP_cell_col_idx)
                        SFP_cell_col_idx_ii = SFP_cell_col_idx(SFP_cell_col_idx_i);
                        SFP_all_new_first(SFP_cell_row(SFP_cell_col_idx_ii), SFP_cell_col(SFP_cell_col_idx_ii)) ...
                            = SFP_color_temp(col_i, 1);
                        SFP_all_new_second(SFP_cell_row(SFP_cell_col_idx_ii), SFP_cell_col(SFP_cell_col_idx_ii)) ...
                            = SFP_color_temp(col_i, 2);
                        SFP_all_new_third(SFP_cell_row(SFP_cell_col_idx_ii), SFP_cell_col(SFP_cell_col_idx_ii)) ...
                            = SFP_color_temp(col_i, 3);
                        
                        % 0.5 means the marked cells;
                        SFP_all_new_label(SFP_cell_row(SFP_cell_col_idx_ii), SFP_cell_col(SFP_cell_col_idx_ii)) = 0.5;
                    end
                end
                
                % the last color;
                SFP_cell_col_idx = find(SFP_cell_col >= SFP_cell_col_min + SFP_cell_col_inter * (SFP_color_num - 1) ...
                    & SFP_cell_col <= SFP_cell_col_max);
                % marked cells;
                for SFP_cell_col_idx_i = 1:length(SFP_cell_col_idx)
                    SFP_cell_col_idx_ii = SFP_cell_col_idx(SFP_cell_col_idx_i);
                    SFP_all_new_first(SFP_cell_row(SFP_cell_col_idx_ii), SFP_cell_col(SFP_cell_col_idx_ii)) = SFP_color_temp(end, 1);
                    SFP_all_new_second(SFP_cell_row(SFP_cell_col_idx_ii), SFP_cell_col(SFP_cell_col_idx_ii)) = SFP_color_temp(end, 2);
                    SFP_all_new_third(SFP_cell_row(SFP_cell_col_idx_ii), SFP_cell_col(SFP_cell_col_idx_ii)) = SFP_color_temp(end, 3);
                    
                    % 0.5 means the marked cells;
                    SFP_all_new_label(SFP_cell_row(SFP_cell_col_idx_ii), SFP_cell_col(SFP_cell_col_idx_ii)) = 0.5;
                end
                
            end
        end
    end
    
    blankMark = sum(SFP_all_new_label == 1,'all');
    SFP_all_new_first(SFP_all_new_label == 1) = color_blankCell(1); % blank cells;
    SFP_all_new_second(SFP_all_new_label == 1) = color_blankCell(2);
    SFP_all_new_third(SFP_all_new_label == 1) = color_blankCell(3);
    
    SFP_all_new_first(SFP_all_new_label == 0) = color_background(1); % 0 means background;
    SFP_all_new_second(SFP_all_new_label == 0) = color_background(2);
    SFP_all_new_third(SFP_all_new_label == 0) = color_background(3);
    
    SFP_all_new(:,:,1) = SFP_all_new_first;
    SFP_all_new(:,:,2) = SFP_all_new_second;
    SFP_all_new(:,:,3) = SFP_all_new_third;
    
    % figure;
    imagesc(SFP_all_new);
    hold on;
    xlim([0 size(SFP_all_new,2)]);
    ylim([0 size(SFP_all_new,1)]);
    
    %% legend;
    if islegend
        for i = 1:length(SFP_mark)
            plotOne(SFP_color{i});
        end
        if blankMark
            plotOne(color_blankCell);
        end
    end
    
    %% contour three; slower; single cell; all cells;
    for i = 1:size(SFP_all,1)
        SFP_all_temp_i = SFP_all(i,:,:);
        SFP_all_temp_i = permute(SFP_all_temp_i,[2 3 1]);
        [~,n,o] = size(SFP_all);
        SFP_all_contour = zeros(n,o);
        %     SFP_mark_contour = zeros(size(SFP_all,2:3));
        SFP_all_contour(SFP_all_temp_i ~= 0) = 1;
        contour(SFP_all_contour,[0.01 0.01],'color',color_cellContour,'LineWidth',1);
    end
    
    if islegend
        if blankMark
            SFP_name = [SFP_name;"other cells"];
        end
        legend(SFP_name,'Interpreter','none');
        legend('boxoff');
    end
    
    hold off;
    axis off;
    
end

function plotOne(PC_color)
    plotmarker = plot(0,0,'.','Color',PC_color,'LineWidth',20);
    set(plotmarker,'visible','off');
end
%% cellDetect.m
% This code is used to detect cells with specific spatial properties,
% including place cell, grid cell, head direction cell, speed cell and may be more.
%
% Input:
%       do_cell_analysis: a struct, choose which cell types to analysis (bool);
%                         (default is 0);
%                         do_cell_analysis.place_cell, 
%                         do_cell_analysis.grid_cell,
%                         do_cell_analysis.border_cell,
%                         do_cell_analysis.head_direction_cell,
%                         do_cell_analysis.speed_cell;
%       cell_name: a struct, custom names of each cell type (char);
%       calcium_time: a m*1 vector, the timestamp of calcium imaging;
%                     (in second, calcium_time(1) = 0 is better);
%       calcium_event: a m*n matrix, n represents the number of cells;
%       pos: a p*3 matrix, [time x y]; (cannot be empty and has three column even in one-dimension)
%       hd_dir: a p*1 vector (degree, 0-360), contains all head direcrions;
%       save_folder: a string, the path that save files;
%       p: a struct, containing necessary paramaters used in cell detection;
%          p.map: a struct, containing necessary paramaters used in map calculation;
%                 will use the default parameters if it is empty;
%                 (details seen in map.m of BNT package)
%          p.field: a struct, containing necessary paramaters used in field calculation;
%                   (details seen in placefield.m or placefield1D.m of BNT package)
%          p.grid: a struct, containing necessary paramaters used in gridness score calculation;
%                  (details seen in gridnessScore.m of BNT package)
%          p.grid2: a struct, containing necessary paramaters used in gridness score calculation;
%                   (properties: 'radii')
%          p.border: a struct, containing necessary paramaters used in border score calculation;
%                    (details seen in borderScore.m of BNT package)
%          p.hd: a struct, containing necessary paramaters used in turning curve calculation;
%                (details seen in turningCurve.m of BNT package)
%          p.hd2: a struct, containing necessary paramaters used in head direction calculation and plot;
%                 (properties: 'sampleTime', 'percentile', 'trajectoryPlot')
%          p.speed: a struct, containing necessary paramaters used in speed score calculation;
%                   (details seen in speedScore.m of BNT package)
%          p.speed2: a struct, containing necessary paramaters used in speed calculation and plot;
%                    (properties:
%                    'speed_firing_rate': a p*n matrix, used for speed cell calculation;
%                    'speed_input': a p*1 vector, used for speed cell calculation;
%                                   it has the same rows with speed_firing_rate;
%                                   can be animal's speed, head direction velocity, etc.
%                                   if speed_firing_rate and speed_input is empty,
%                                   the script will process default calculation;
%                    'speed_sample_time': the sample time of speed_firing_rate and speed_input;
%                                         (default: 0.02 (sec));
%                    'speed_win': the smooth window of speed calculation;
%                                 (default: 100 (frame));
%                    'speed_range': a 1*3 vector, [do_speed_range, min_speed, max_speed];
%                                   the range of speed, the value below it or above it will be set to nan;
%                                   (default: [1 2 0]);
%                    'speed_bin': the bin of speed tuning plots;
%                                 (default: 2 (cm/s));)
%          p.behav_limit: a 1*4 vector, [xmin xmax ymin ymax], will be used in trajectory plots;
%                         (same as p.map.limits)
%
% Optional input:
%       min_event_num: the minimum number of events;
%       spatial_dimension: a scale or a 1*2 vector;
%                          the first digit is 1 or 2, representing linear or two-dimension;
%                          the second digit is the column used to measure or left in empty;
%                          e.g. [1 1] or [1 2] or 2, (default: 2);
%       speed_filter: a 1*3 vector, [do_speed_filter, min_speed_filter, max_speed_filter];
%                     (default: [1 2.5 0]);
%       do_shuffle: bool, do shuffle process or not (default: 1);
%       shuffle_num: the number of perform shuffles (default: 1000);
%       threshold_prc: a number, used to determine the threshold of shuffled data;
%                      (from 1-99, default: 95);
%                      OR a struct, containing the threshold of each cell type;
%       draw_fig: bool, whether draw plots or not (default: 1);
%                 OR a struct, containing whrther each types;
%       SFP: a cellNum*width*height matrix, spatial footprint;
%       fig_fmt: a 1*n cell, containing all the formats of saved figures;
%                (default: {'.fig'});
%       color_cell: a struct, containing all the colors of cells in SFP figures ([R G B]);
%                   (default is [0.7 0 0]);
%                   color_cell.cell;
%                   color_cell.place_cell, 
%                   color_cell.grid_cell,
%                   color_cell.border_cell,
%                   color_cell.head_direction_cell,
%                   color_cell.speed_cell;
%
% Output:
%
% Usage: details seen in cellDetect_run.m.
%
%==========================================================================
% v1.0: March, 2023.
% Users can choose which cell type to calculate and detect.
% The rate map and trajectories of place cell, grid cell and border cell
% are plotted in the same folder.
% The calculation and plots of speed cell are not completed.
%
% v1.1: April, 2023.
% Users can choose the prctile of shuffle for each cell type;
% Add head direction stability.
%
% v2.0: May, 2023.
% Choose to do shuffle process or not;
% Choose which figure to draw or not;
% Change the struct of inputs, cell detect and save results parts;
% Add plotting cell path with head direction function.
%
% v2.2: June, 2023.
% Correct the part of detecting speed cell with the help of Emilio.
% The parameters and inputs must be treated carefully.
%
% v2.3: July, 2023.
% Add linear track calculation (one dimension, v2.1 not tested);
% Change the order of head_direction_cell and border_cell;
% Modify the calculation in function cellStability;
% Integrate the input of speed calculation into p.speed2;
% Save parameters used in the calculation;
% Add speed tuning curve in the output.
%
% v3.0: August, 2023.
% Add the min_event_num parameter;
% Test and fix bugs.
%
% v3.1: November, 2023.
% Fix bugs of detecting speed cell.
%
% Created by Xiang Zhang, March, 2023.

function cellDetect(do_cell_analysis, cell_name, calcium_time, calcium_event, pos, hd_dir, save_folder, p, varargin)
    %% check inputs and optional inputs;
    do_cell_analysis = doCellAnalysis(do_cell_analysis);
    if isempty(do_cell_analysis), disp('No cell type need to be detected, please check the inputs.'); return; end
    
    cell_name = cellName(cell_name);
    
    if isempty(pos), error('No valid data in position, please check it.'); end
    if size(pos, 2) ~= 3, error('The position data is not correct in size, please check it.'); end
    if do_cell_analysis.head_direction_cell && isempty(hd_dir), ...
            error('No valid data in head direction, please check it.'); end
    if do_cell_analysis.head_direction_cell && min(hd_dir) < 0, hd_dir = mod(hd_dir, 360); end
    
    inp = inputParser;
    addRequired(inp,'do_cell_analysis', @isstruct);
    addRequired(inp, 'cell_name', @isstruct);
    addRequired(inp, 'calcium_time', @isvector);
    addRequired(inp, 'calcium_event', @ismatrix);
    addRequired(inp, 'pos', @ismatrix);
    addRequired(inp, 'hd_dir'); % , @isvector
    addRequired(inp, 'save_folder', @ischar);
    addRequired(inp, 'p', @isstruct);
    
    addParameter(inp, 'min_event_num', 10);
    addParameter(inp, 'spatial_dimension', 2);
    addParameter(inp, 'speed_filter', [1, 2.5, 0]);
    addParameter(inp, 'do_shuffle', 1);
    addParameter(inp, 'shuffle_num', 1000);
    addParameter(inp, 'threshold_prc', 95);
    addParameter(inp, 'draw_fig', 1);
    addParameter(inp, 'SFP', []);
    addParameter(inp, 'fig_fmt', {'.fig'});
    addParameter(inp, 'color_cell', struct, @isstruct);
    
    parse(inp, do_cell_analysis, cell_name, calcium_time, calcium_event, pos, hd_dir, save_folder, p, varargin{:});
    min_event_num = inp.Results.min_event_num;
    spatial_dimension = inp.Results.spatial_dimension;
    speed_filter = inp.Results.speed_filter; if speed_filter(3) == 0, speed_filter(3) = Inf; end
    do_shuffle = inp.Results.do_shuffle;
    shuffle_num = inp.Results.shuffle_num;
    threshold_prc = inp.Results.threshold_prc;
    draw_fig = inp.Results.draw_fig;
    SFP = inp.Results.SFP;
    fig_fmt = inp.Results.fig_fmt;
    color_cell = inp.Results.color_cell;
    
    % sptial dimension;
    if spatial_dimension(1) == 1
        if length(spatial_dimension) ~= 2, error('The input ''spatial_dimension'' is not valid, please check it.'); end
        do_cell_analysis.grid_cell = 0;
        do_cell_analysis.border_cell = 0;
        % do_cell_analysis.head_direction_cell = 0;
    end
    
    % figure;
    draw_fig = drawFig(draw_fig, do_cell_analysis);
    if draw_fig.cell && isempty(SFP), error('No valid data in SFP, please check it.'); end
    if draw_fig.cell, color_cell = selectColor(color_cell); end
    if (draw_fig.cellpath || draw_fig.merge) && ~isfield(p, 'behav_limit')
        error('No valid data in behav_limit, please check it.'); end
    
    % threshold;
    if shuffle_num == 0, do_shuffle = 0; end
    if do_shuffle == 0, shuffle_num = 0; end
    if do_shuffle, threshold_prc = thresholdPrctile(threshold_prc); end
    
    %% initialization;
    % make folders;
    makeFolders(do_cell_analysis, draw_fig, cell_name, save_folder);
    
    % parameters;
    if ~isfield(p, 'map'), p.map = struct; end, para.map = collectPara(p.map);
    if ~isfield(p, 'field'), p.field = struct; end, para.field = collectPara(p.field);
    para = pInit(do_cell_analysis, p, para);
    numNeurons = size(calcium_event,2);
    
    % cells;
    [cell_all, output_shuffle] = shuffleInit(do_cell_analysis, numNeurons, shuffle_num);
    
    % speed filter;
    [pos_modified, pos_modified_full, hd_dir_modified] = speedFilter(speed_filter, spatial_dimension, pos, hd_dir);
    
    % stability;
    splitDataArray_modified = splitData(pos_modified);
    splitDataArray_hd_modified = splitData([pos(:,1), hd_dir_modified]);
    
    % head direction cell;
    if do_cell_analysis.head_direction_cell
        trajectory_map = general.circHist(hd_dir, p.hd.binWidth);
        trajectory_map = general.circSmooth(trajectory_map, p.hd.smooth);
        p.hd2.trajectory_map = trajectory_map;
        p.hd2.trajectory_hd = trajectory_map / max(trajectory_map);
    end
    
    % speed cell;
    if do_cell_analysis.speed_cell
        p = speedInit(do_cell_analysis, p);
        p.speed2 = speedProcess(do_cell_analysis, calcium_event, calcium_time, pos(:,1), pos, p.speed2);
    end
    
    %% main function;
    % if exist([save_folder, '/output.mat'], 'file'), load([save_folder, '/output.mat'], 'output');
    % else
    output = table('Size',[numNeurons 18], ...
        'VariableTypes',["cell", "cell", "cell", "cell", "double", "double", ...
        "cell", "double", "cell", "double", "double", "double", ...
        "cell", "cell", "double", "cell", "cell", "double"], ...
        'VariableNames',["map", "fields_map", "fields", "map_stability", "information_content", "information_rate", ...
        "autocorrelogram", "grid_score", "grid_stat", "center_field", "score_radius", "border_score", ...
        "tc", "tcStat", "curMVL", "angle_stability", "speed_tuning_curve", "speed_score"]);
    % end
    
    parfor cell_i = 1:numNeurons
        warning off all;
        k = calcium_event(:,cell_i) > 0; % 3 * std(DeconvSignals(:,i))
        if sum(k) < min_event_num, continue; end
        
        spike_pos = calcium_time(k); %#ok<PFBNS> % the time of spike
        spike_pos = spikePos(spike_pos, pos_modified(:,1), pos_modified_full, hd_dir_modified); %#ok<PFBNS>
        if size(spike_pos,1) < min_event_num, continue; end % isempty(spike_pos)
        
        %% basic information;
        output_temp = cellInfo(do_cell_analysis, cell_i, spike_pos, splitDataArray_modified, splitDataArray_hd_modified, ...
            pos_modified, hd_dir_modified, spatial_dimension, p, para);
        output(cell_i,:) = output_temp; %#ok<*PFOUS>
        
        %% shuffle;
        output_shuffle_temp = shuffleProcess(do_shuffle, do_cell_analysis, cell_i, k, calcium_time, ...
            shuffle_num, threshold_prc, output_temp, pos_modified_full, hd_dir_modified, spatial_dimension, ...
            splitDataArray_modified, splitDataArray_hd_modified, p, para);
        
        %% cell detect and draw plots;
        [cell_all(cell_i,:), output_shuffle(cell_i,:), cell_sig] = ...
            cellDetectFcn(do_cell_analysis, cell_i, output_temp, output_shuffle_temp);
        cellPropPlot(draw_fig, do_cell_analysis, spatial_dimension, cell_i, output_temp, SFP, pos, spike_pos, ...
            cell_sig, cell_name, color_cell, p, save_folder, fig_fmt);
        
    end
    
    %% information collect and storage;
    save([save_folder, '/parameters.mat'], ...
        'do_cell_analysis', 'cell_name', 'calcium_time', 'calcium_event', 'pos', 'hd_dir', 'save_folder', 'spatial_dimension', ...
        'speed_filter', 'do_shuffle', 'shuffle_num', 'threshold_prc', 'draw_fig', 'SFP', 'fig_fmt', 'color_cell');
    save([save_folder, '/p.mat'], 'p');
    save([save_folder, '/output.mat'], 'output', '-v7.3');
    % save([save_folder, '/output_shuffle.mat'], 'output_shuffle', '-v7.3');
    
    saveCellResults(do_cell_analysis, do_shuffle, output, output_shuffle, cell_all, cell_name, ...
        save_folder, draw_fig, SFP, fig_fmt, color_cell, p);
    
    %% end;
    close all;
    disp([save_folder, ' finished calculation and files were saved.']);
    
end

%% functions;
%-------------------------------------------------------------------------%
%                              parameters                                 %
%-------------------------------------------------------------------------%
% determine cell types;
function do_cell_analysis_output = doCellAnalysis(do_cell_analysis)
    do_cell_analysis_output = struct;
    do_cell_analysis_output.place_cell = 0;
    do_cell_analysis_output.grid_cell = 0;
    do_cell_analysis_output.border_cell = 0;
    do_cell_analysis_output.head_direction_cell = 0;
    do_cell_analysis_output.speed_cell = 0;
    
    cell_types = fieldnames(do_cell_analysis);
    return_sum = 0;
    if isempty(cell_types)
    else
        for type_i = 1:length(cell_types)
            do_cell_analysis_output.(cell_types{type_i}) = do_cell_analysis.(cell_types{type_i});
%             switch cell_types{i}
%                 case 'place_cell'
%                     do_cell_analysis_output.place_cell = do_cell_analysis.(cell_types{i});
%                 case 'grid_cell'
%                     do_cell_analysis_output.grid_cell = do_cell_analysis.(cell_types{i});
%                 case 'border_cell'
%                     do_cell_analysis_output.border_cell = do_cell_analysis.(cell_types{i});
%                 case 'head_direction_cell'
%                     do_cell_analysis_output.head_direction_cell = do_cell_analysis.(cell_types{i});
%                 case 'speed_cell'
%                     do_cell_analysis_output.speed_cell = do_cell_analysis.(cell_types{i});
%             end
            return_sum = return_sum + do_cell_analysis.(cell_types{type_i});
        end
    end
    if ~return_sum, do_cell_analysis_output = []; end
end

function cell_name_output = cellName(cell_name)
    cell_name_output = struct;
    cell_name_output.place_cell = 'place_cell';
    cell_name_output.grid_cell = 'grid_cell';
    cell_name_output.border_cell = 'border_cell';
    cell_name_output.head_direction_cell = 'head_direction_cell';
    cell_name_output.speed_cell = 'speed_cell';
    
    cell_types = fieldnames(cell_name);
    if ~isempty(cell_types)
        for type_i = 1:length(cell_types)
            if ~isempty(cell_name.(cell_types{type_i}))
                cell_name_output.(cell_types{type_i}) = cell_name.(cell_types{type_i});
            end
        end
    end
end

function color_cell_output = selectColor(color_cell)
    color_cell_output = struct;
    color_cell_output.cell = [0.8 0 0];
    color_cell_output.place_cell = [0.8 0 0];
    color_cell_output.grid_cell = [0.8 0 0];
    color_cell_output.border_cell = [0.8 0 0];
    color_cell_output.head_direction_cell = [0.8 0 0];
    color_cell_output.speed_cell = [0.8 0 0];
    
    cell_types = fieldnames(color_cell);
    if ~isempty(cell_types)
        for type_i = 1:length(cell_types)
            if ~isempty(color_cell.(cell_types{type_i}))
                color_cell_output.(cell_types{type_i}) = color_cell.(cell_types{type_i});
            end
        end
    end
end

function draw_fig_output = drawFig(draw_fig, do_cell_analysis)
    [draw_fig_output.cell, draw_fig_output.ratemap, ...
        draw_fig_output.cellpath, draw_fig_output.autocorrelogram, ...
        draw_fig_output.turningcurve, draw_fig_output.cellpath_hd, ...
        draw_fig_output.merge, draw_fig_output.speedtuning] = deal(1);
    if isnumeric(draw_fig)
        [draw_fig_output.cell, draw_fig_output.ratemap, ...
            draw_fig_output.cellpath, draw_fig_output.autocorrelogram, ...
            draw_fig_output.turningcurve, draw_fig_output.cellpath_hd, ...
            draw_fig_output.merge, draw_fig_output.speedtuning] = deal(draw_fig);
    elseif isstruct(draw_fig)
        draw_types = fieldnames(draw_fig);
        for type_i = 1:length(draw_types)
            switch draw_types{type_i}
                case {'cell', 'ratemap', 'cellpath', 'autocorrelogram', ...
                        'turningcurve', 'cellpath_hd', 'merge', 'speedtuning'}
                    draw_fig_output.(draw_types{type_i}) = draw_fig.(draw_types{type_i});
                otherwise
                    error('Please cheack the input ''draw_fig''');
            end
        end
    end
    
    if ~do_cell_analysis.place_cell && ~do_cell_analysis.grid_cell && ...
            ~do_cell_analysis.border_cell && ~do_cell_analysis.head_direction_cell
        draw_types = {'ratemap', 'cellpath', 'autocorrelogram', 'turningcurve', 'cellpath_hd', 'merge'};
        for type_i = 1:length(draw_types)
            draw_fig_output.(draw_types{type_i}) = 0;
        end
    end
end

function threshold_prc_output = thresholdPrctile(threshold_prc)
    [threshold_prc_output.map_stability_half, threshold_prc_output.map_stability_time, ...
        threshold_prc_output.information_content, threshold_prc_output.information_rate, ...
        threshold_prc_output.grid_score, threshold_prc_output.border_score, threshold_prc_output.curMVL, ...
        threshold_prc_output.angle_stability_half, threshold_prc_output.angle_stability_time, ...
        threshold_prc_output.speed_score] = deal(95);
    if isnumeric(threshold_prc)
        [threshold_prc_output.map_stability_half, threshold_prc_output.map_stability_time, ...
            threshold_prc_output.information_content, threshold_prc_output.information_rate, ...
            threshold_prc_output.grid_score, threshold_prc_output.border_score, threshold_prc_output.curMVL, ...
            threshold_prc_output.angle_stability_half, threshold_prc_output.angle_stability_time, ...
            threshold_prc_output.speed_score] = deal(threshold_prc);
    elseif isstruct(threshold_prc)
        shuffle_types = fieldnames(threshold_prc);
        for type_i = 1:length(shuffle_types)
            switch shuffle_types{type_i}
                case 'place_cell'
                    [threshold_prc_output.map_stability_half, threshold_prc_output.map_stability_time, ...
                        threshold_prc_output.information_content, threshold_prc_output.information_rate] = ...
                        deal(threshold_prc.(shuffle_types{type_i}));
                case 'grid_cell'
                    [threshold_prc_output.map_stability_half, threshold_prc_output.map_stability_time, ...
                        threshold_prc_output.grid_score] = deal(threshold_prc.(shuffle_types{type_i}));
                case 'border_cell'
                    [threshold_prc_output.map_stability_half, threshold_prc_output.map_stability_time, ...
                        threshold_prc_output.border_score] = deal(threshold_prc.(shuffle_types{type_i}));
                case 'head_direction_cell'
                    [threshold_prc_output.curMVL, threshold_prc_output.angle_stability_half, ...
                        threshold_prc_output.angle_stability_time] = deal(threshold_prc.(shuffle_types{type_i}));
                case 'speed_cell'
                    [threshold_prc_output.speed_score] = deal(threshold_prc.(shuffle_types{type_i}));
                case {'map_stability_half', 'map_stability_time', ...
                        'information_content', 'information_rate', 'grid_score', 'border_score', ...
                        'curMVL', 'angle_stability_half', 'angle_stability_time', 'speed_score'}
                    threshold_prc_output.(shuffle_types{type_i}) = threshold_prc.(shuffle_types{type_i});
                otherwise
                    error('Please cheack the input ''threshold_prc''');
            end
        end
    end
end

function makeFolders(do_cell_analysis, draw_fig, cell_name, save_folder) % need modifications;
    makeFolder(save_folder)
    cell_types = fieldnames(do_cell_analysis);
    for type_i = 1:length(cell_types)
        if do_cell_analysis.(cell_types{type_i})
            makeFolder(strcat(save_folder, '/', cell_name.(cell_types{type_i}))); % shuffles
            switch cell_types{type_i}
                case 'place_cell'
                    
                case 'grid_cell'
                    if draw_fig.autocorrelogram
                        makeFolder(strcat(save_folder, '/', cell_name.(cell_types{type_i}), '/autocorrelogram')); end
                case 'border_cell'
                    
                case 'head_direction_cell'
                    if draw_fig.turningcurve
                        makeFolder(strcat(save_folder, '/', cell_name.(cell_types{type_i}), '/turningcurve')); end
                    if draw_fig.cellpath_hd
                        makeFolder(strcat(save_folder, '/', cell_name.(cell_types{type_i}), '/cellpath')); end
                case 'speed_cell'
                    if draw_fig.speedtuning
                        makeFolder(strcat(save_folder, '/', cell_name.(cell_types{type_i}), '/speedtuning')); end
            end
        end
    end
    if do_cell_analysis.place_cell || do_cell_analysis.grid_cell || do_cell_analysis.border_cell
        makeFolder(strcat(save_folder, '/map'));
        if draw_fig.ratemap, makeFolder(strcat(save_folder, '/ratemap')); end
        if draw_fig.cellpath, makeFolder(strcat(save_folder, '/cellpath')); end
        if draw_fig.merge, makeFolder(strcat(save_folder, '/merge')); end
    end
    % mkdir(strcat(save_folder, '/mapstability'));
    if draw_fig.cell, makeFolder(strcat(save_folder, '/cell')); end
end

function makeFolder(folder_name)
    if ~exist(folder_name, 'dir'), mkdir(folder_name); end
end

function para = collectPara(p)
    para = {};
    if isempty(p), return; end
    p_var = fieldnames(p);
    for var_i = 1:length(p_var)
        para{1,end+1} = p_var{var_i}; %#ok<AGROW>
        para{1,end+1} = p.(p_var{var_i}); %#ok<AGROW>
    end
end

function para = pInit(do_cell_analysis, p, para)
    if do_cell_analysis.grid_cell
        if ~isfield(p, 'grid'), p.grid = struct; end, para.grid = collectPara(p.grid);
        if ~isfield(p, 'grid2'), p.grid2.radii = [5 5]; else, if ~isfield(p.grid2, 'radii'), p.grid2.radii = [5 5]; end, end
    end
    
    if do_cell_analysis.border_cell
        if ~isfield(p, 'border'), p.border = struct; end, para.border = collectPara(p.border);
    end
    
    if do_cell_analysis.head_direction_cell
        if ~isfield(p, 'hd'), p.hd = struct; end, para.hd = collectPara(p.hd);
    end
    
    if do_cell_analysis.speed_cell
        if ~isfield(p, 'speed'), p.speed = struct; end, para.speed = collectPara(p.speed);
    end
end

function [cell_all, output_shuffle] = shuffleInit(do_cell_analysis, numNeurons, shuffle_num)
    cell_all = table('Size',[numNeurons, 5], ...
        'VariableTypes',["double", "double", "double", "double", "double"], ...
        'VariableNames',["place_cell", "grid_cell", "border_cell", "head_direction_cell", "speed_cell"]);
    output_shuffle = table('Size',[numNeurons, 10], ...
        'VariableTypes',["double", "double", "double", "double", ...
        "double", "double", "double", "double", "double", "double"], ...
        'VariableNames',["map_stability_half_shuffle", "map_stability_time_shuffle", ...
        "information_content_shuffle", "information_rate_shuffle", "grid_score_shuffle", ...
        "border_score_shuffle", "curMVL_shuffle", ...
        "angle_stability_half_shuffle", "angle_stability_time_shuffle", "speed_score_shuffle"]);
    
    if do_cell_analysis.place_cell || do_cell_analysis.grid_cell || do_cell_analysis.border_cell
        [output_shuffle.map_stability_half_shuffle, output_shuffle.map_stability_time_shuffle] = ...
            deal(nan(numNeurons, shuffle_num + 3));
    end
    
    if do_cell_analysis.place_cell
        % information_all = table(nan(numNeurons,1), nan(numNeurons,1), 'VariableNames',["content", "rate"]);
        [output_shuffle.information_content_shuffle, output_shuffle.information_rate_shuffle] = ...
            deal(nan(numNeurons, shuffle_num + 3));
        cell_all.place_cell = nan(numNeurons,1);
    end
    
    if do_cell_analysis.grid_cell
        output_shuffle.grid_score_shuffle = nan(numNeurons, shuffle_num + 3);
        cell_all.grid_cell = nan(numNeurons,1);
    end
    
    if do_cell_analysis.border_cell
        output_shuffle.border_score_shuffle = nan(numNeurons, shuffle_num + 3);
        cell_all.border_cell = nan(numNeurons,1);
    end
    
    if do_cell_analysis.head_direction_cell
        [output_shuffle.curMVL_shuffle, ...
            output_shuffle.angle_stability_half_shuffle, output_shuffle.angle_stability_time_shuffle] = ...
            deal(nan(numNeurons, shuffle_num + 3));
        cell_all.head_direction_cell = nan(numNeurons,1);
    end
    
    if do_cell_analysis.speed_cell
        output_shuffle.speed_score_shuffle = nan(numNeurons, shuffle_num + 3);
        cell_all.speed_cell = nan(numNeurons,1);
    end
end

function p = speedInit(do_cell_analysis, p)
    if ~do_cell_analysis.speed_cell, return; end
    p_input = p;
    % defalut values;
    p.speed2.speed_firing_rate = [];
    p.speed2.speed_input = [];
    p.speed2.speed_sample_time = 0.2;
    p.speed2.speed_win = 100;
    p.speed2.speed_range = [1, 2, 0];
    p.speed2.speed_bin = 2;
    
    % input values;
    speed_input = fieldnames(p_input.speed2);
    if isempty(speed_input), return; end
    for input_i = 1:length(speed_input)
        p.speed2.(speed_input{input_i}) = p_input.speed2.(speed_input{input_i});
    end
    if p.speed2.speed_range(3) == 0 || isinf(p.speed2.speed_range(3))
        p.speed2.speed_range(3) = Inf;
        speed_range = p.speed2.speed_range(2):p.speed2.speed_bin:max(p_input.speed2.speed_input);
        p.speed2.speed_range_bin = [speed_range, speed_range(end)+p.speed2.speed_bin];
    else
        speed_range = p.speed2.speed_range(2):p.speed2.speed_bin:p.speed2.speed_range(3);
        p.speed2.speed_range_bin = [speed_range, speed_range(end)+p.speed2.speed_bin];
    end
end


%-------------------------------------------------------------------------%
%                             data process                                %
%-------------------------------------------------------------------------%
function [pos_modified, pos_modified_full, hd_dir_modified] = speedFilter(speed_filter, spatial_dimension, pos, hd_dir)
    pos_modified_full = pos; % includes t, x, y;
    hd_dir_modified = hd_dir;
    % angular head velocity;
%     if do_cell_analysis.angular_head_velocity
%         ahv_modified = analyses.angularHeadVelocity([pos(:,1), hd_dir]);
%     else, ahv_modified = [];
%     end
    
    if ~speed_filter(1)
        if spatial_dimension(1) == 2, pos_modified = pos;
        elseif spatial_dimension(1) == 1, pos_modified = pos(:, [1 spatial_dimension(2)+1]); end % select x or y;
    else
        speed = speed2D(pos(:,2), pos(:,3), pos(:,1));
        speed_filter_i = speed < speed_filter(2) | speed > speed_filter(3);
        pos_modified_full(speed_filter_i, 2:end) = nan;
        if spatial_dimension(1) == 2, pos_modified = pos_modified_full;
        elseif spatial_dimension(1) == 1, pos_modified = pos_modified_full(:, [1 spatial_dimension(2)+1]); end % select x or y;
        if ~isempty(hd_dir), hd_dir_modified(speed_filter_i, :) = nan; end
%         if do_cell_analysis.angular_head_velocity, ahv_modified(speed_filter_i, :) = nan; end
    end
end

% split data in two ways;
function [splitDataArray, pos_first, pos_second, pos_odd, pos_even] = splitData(pos)
    % spilt data in four ways;
    % 1 Row:    First half of data (Half and half type split)
    % 2 Row:    Second half of data (Half and half type split)
    % 3 Row:    First half of data (Binned type split) (1 min and 1 min)
    % 4 Row:    Second half of data (Binned type split) (1 min and 1 min)
    splitDataArray = cell(4,1);
    
    duration = pos(end,1); % msec; - pos(1,1);
    half_ind = pos(:, 1) <= duration/2; % + pos(1,1));
    [splitDataArray{1,1}, pos_first] = deal(pos(half_ind, :));
    
    half_ind = pos(:, 1) > duration/2; % + pos(1,1));
    [splitDataArray{2,1}, pos_second] = deal(pos(half_ind, :));
    
    splitDataArray(3:4, :) = splitOddEven(pos);
    pos_odd = splitDataArray{3,1};
    pos_even = splitDataArray{4,1};
end

% split data in odd and even times;
function splitDataArray = splitOddEven(pos)
    duration = pos(end,1); %  - pos(1,1);
    [txy3, txy4] = deal(zeros(size(pos,1), size(pos,2)));
    [totSamp3,totSamp4] = deal(0);
    
    nBins = ceil(duration / 60);
    t_start = 0; % pos(1,1);
    t_stop = t_start + 60;
    for spilt_i = 1:nBins
        ind = find(pos(:,1) >= t_start & pos(:,1) < t_stop);
        samps = length(ind);
        if samps > 0
            if mod(spilt_i,2) % odd;
                txy3(totSamp3+1:totSamp3+samps, :) = pos(ind, :);
                totSamp3 = totSamp3 + samps;
            else % even;
                txy4(totSamp4+1:totSamp4+samps, :) = pos(ind, :);
                totSamp4 = totSamp4 + samps;
            end
        end
        t_start = t_stop;
        t_stop = t_start + 60;
    end
    
    splitDataArray{1,1} = txy3(1:totSamp3, :);
    splitDataArray{2,1} = txy4(1:totSamp4, :);
end

function spike_pos = spikePos(spike_pos, behav_time, pos, hd_dir)
    if isempty(spike_pos) || (isempty(pos) && isempty(hd_dir)), return; end
    
    spkInd = knnsearch(behav_time, spike_pos(:,1));
    if ~isempty(pos), spike_pos(:, 2:3) = pos(spkInd, 2:3); else, spike_pos(:, 2:3) = nan; end
    if ~isempty(hd_dir), spike_pos(:, 4) = hd_dir(spkInd); end
    
    spike_pos(isnan(spike_pos(:,2)), 3) = nan;
    spike_pos(isnan(spike_pos(:,3)), 2) = nan;
    spike_pos(sum(isnan(spike_pos(:,2:end)),2) == (size(spike_pos,2) - 1), :) = [];
%     spike_pos(isnan(spike_pos(:,2)),:) = []; % delete before use it;
end

function map_stability = cellStability(spike_pos, splitDataArray_modified, para_map)
    map_stability = struct;
    map_stability = cellStability_half(map_stability, spike_pos, splitDataArray_modified, para_map);
    map_stability = cellStability_time(map_stability, spike_pos, splitDataArray_modified, para_map);
end

% half stability;
function map_stability = cellStability_half(map_stability, spike_pos, splitDataArray_modified, para_map)
    spike_pos_first = spike_pos(spike_pos(:,1) <= splitDataArray_modified{1}(end,1), :);
    spike_pos_second = spike_pos(spike_pos(:,1) > splitDataArray_modified{1}(end,1), :);
    
    pos_first_modified = splitDataArray_modified{1};
    pos_second_modified = splitDataArray_modified{2};
    
    if ~isempty(spike_pos_first) && ~isempty(spike_pos_second)
        map_first = analyses.map(pos_first_modified, spike_pos_first(:,1), para_map{:});
        map_first.spike_pos = spike_pos_first;
        map_second = analyses.map(pos_second_modified, spike_pos_second(:,1), para_map{:});
        map_second.spike_pos = spike_pos_second;
        
        % correlation;
        r_half = zeroLagCorrelation(map_first.z,map_second.z);
    else
        map_first = [];
        map_second = [];
        r_half = nan;
    end
    map_stability.map_first = map_first;
    map_stability.map_second = map_second;
    map_stability.half_correlation = r_half;
end

% time binned stability;
function map_stability = cellStability_time(map_stability, spike_pos, splitDataArray_modified, para_map)
    
    splitDataArray_spike = splitOddEven(spike_pos);
    
    spike_pos_odd = splitDataArray_spike{1};
    spike_pos_even = splitDataArray_spike{2};
    
    if ~isempty(spike_pos_odd) && ~isempty(spike_pos_even)
        pos_odd_modified = splitDataArray_modified{3};
        pos_even_modified = splitDataArray_modified{4};
        
        map_odd = analyses.map(pos_odd_modified, spike_pos_odd(:,1), para_map{:});
        map_odd.spike_pos = spike_pos_odd;
        map_even = analyses.map(pos_even_modified, spike_pos_even(:,1), para_map{:});
        map_even.spike_pos = spike_pos_even;
        
        % correlation;
        r_time = zeroLagCorrelation(map_odd.z,map_even.z);
    else
        map_odd = [];
        map_even = [];
        r_time = nan;
    end
    map_stability.map_odd = map_odd;
    map_stability.map_even = map_even;
    map_stability.time_correlation = r_time;
end

function angle_stability = cellStability_hd(spike_pos, splitDataArray_hd_modified, p, para)
    angle_stability = struct;
    angle_stability = cellStability_hd_half(angle_stability, spike_pos, splitDataArray_hd_modified, p, para);
    angle_stability = cellStability_hd_time(angle_stability, spike_pos, splitDataArray_hd_modified, p, para);
    
    angle_stability.p.hd = p.hd;
    angle_stability.p.hd2 = p.hd2;
end

% half stability;
function angle_stability = cellStability_hd_half(angle_stability, spike_pos, splitDataArray_hd_modified, p, para)
    spike_pos_first = spike_pos(spike_pos(:,1) <= splitDataArray_hd_modified{1}(end,1), :);
    spike_pos_second = spike_pos(spike_pos(:,1) > splitDataArray_hd_modified{1}(end,1), :);
    spkDir_first = spike_pos_first(:, 2); spkDir_first(isnan(spkDir_first)) = [];
    spkDir_second = spike_pos_second(:, 2); spkDir_second(isnan(spkDir_second)) = [];
    
    hd_dir_first_modified = splitDataArray_hd_modified{1};
    hd_dir_second_modified = splitDataArray_hd_modified{2};
    
    if ~isempty(spkDir_first) && ~isempty(spkDir_second) % && ~isempty(hd_dir_first_modified) && ~isempty(hd_dir_second_modified)
        tc_first = analyses.turningCurve(spkDir_first, hd_dir_first_modified(:,2), p.hd2.sampleTime, para.hd{:});
        tcStat_first = analyses.tcStatistics(tc_first, p.hd.binWidth, p.hd2.percentile);
        tcStat_first.spkDir = spkDir_first;
        
        tc_second = analyses.turningCurve(spkDir_second, hd_dir_second_modified(:,2), p.hd2.sampleTime, para.hd{:});
        tcStat_second = analyses.tcStatistics(tc_second, p.hd.binWidth, p.hd2.percentile);
        tcStat_second.spkDir = spkDir_second;
        
        % correlation;
        r_half = corr(tc_first(:,2), tc_second(:,2));
    else
        tc_first = [];
        tcStat_first = [];
        tc_second = [];
        tcStat_second = [];
        r_half = nan;
    end
    angle_stability.tc_first = tc_first;
    angle_stability.tcStat_first = tcStat_first;
    angle_stability.tc_second = tc_second;
    angle_stability.tcStat_second = tcStat_second;
    angle_stability.half_correlation = r_half;
end

% time binned stability;
function angle_stability = cellStability_hd_time(angle_stability, spike_pos, splitDataArray_hd_modified, p, para)
    
    splitDataArray_spike = splitOddEven(spike_pos);
    
    spike_pos_odd = splitDataArray_spike{1};
    spike_pos_even = splitDataArray_spike{2};
    spkDir_odd = spike_pos_odd(:, 2); spkDir_odd(isnan(spkDir_odd)) = [];
    spkDir_even = spike_pos_even(:, 2); spkDir_even(isnan(spkDir_even)) = [];
    
    if ~isempty(spkDir_odd) && ~isempty(spkDir_even)
        hd_dir_odd_modified = splitDataArray_hd_modified{3};
        hd_dir_even_modified = splitDataArray_hd_modified{4};
        
        tc_odd = analyses.turningCurve(spkDir_odd, hd_dir_odd_modified(:,2), p.hd2.sampleTime, para.hd{:});
        tcStat_odd = analyses.tcStatistics(tc_odd, p.hd.binWidth, p.hd2.percentile);
        tcStat_odd.spkDir = spkDir_odd;
        
        tc_even = analyses.turningCurve(spkDir_even, hd_dir_even_modified(:,2), p.hd2.sampleTime, para.hd{:});
        tcStat_even = analyses.tcStatistics(tc_even, p.hd.binWidth, p.hd2.percentile);
        tcStat_even.spkDir = spkDir_even;
        
        % correlation;
        r_time = corr(tc_odd(:,2), tc_even(:,2));
    else
        tc_odd = [];
        tcStat_odd = [];
        tc_even = [];
        tcStat_even = [];
        r_time = nan;
    end
    angle_stability.tc_odd = tc_odd;
    angle_stability.tcStat_odd = tcStat_odd;
    angle_stability.tc_even = tc_even;
    angle_stability.tcStat_even = tcStat_even;
    angle_stability.time_correlation = r_time;
end

function p_input = speedProcess(do_cell_analysis, calcium_event, calcium_time, behav_time, pos, p_input)
    if ~do_cell_analysis.speed_cell, p_input.speed_firing_rate = nan(size(calcium_event)); return; end
    
    % generate the same time line for calcium event and speed;
    time_all = 0:p_input.speed_sample_time:min(calcium_time(end), behav_time(end));
    
    % smooth the firing rate;
    if isempty(p_input.speed_firing_rate)
        p_input.speed_firing_rate = interp1(calcium_time, calcium_event, time_all);
        p_input.speed_firing_rate = cell2mat(arrayfun(@(x) ...
            general.smoothGauss(p_input.speed_firing_rate(:,x), 0.4/p_input.speed_sample_time), ... % 20
            1:size(p_input.speed_firing_rate,2), 'UniformOutput', false));
    end
    
    % calculate and smooth speed;
    if isempty(p_input.speed_input)
        % speed_win = 100; % time window;
        pos_smooth = interp1(behav_time, pos, time_all);
        speed_x = pos_smooth(p_input.speed_win+1:end,2) - pos_smooth(1:end-p_input.speed_win,2);
        speed_y = pos_smooth(p_input.speed_win+1:end,3) - pos_smooth(1:end-p_input.speed_win,3);
        p_input.speed_input = [nan(p_input.speed_win/2,1); ...
            sqrt(speed_x.^2 + speed_y.^2) / (p_input.speed_win * p_input.speed_sample_time); nan(p_input.speed_win/2,1)];
    end
    if p_input.speed_range(1)
        p_input.speed_input(p_input.speed_input < p_input.speed_range(2)) = nan;
        p_input.speed_input(p_input.speed_input > p_input.speed_range(3)) = nan;
    end
end

%-------------------------------------------------------------------------%
%                            main calculation                             %
%-------------------------------------------------------------------------%
function output = cellInfo(do_cell_analysis, cell_i, spike_pos, splitDataArray_modified, splitDataArray_hd_modified, ...
        pos_modified, hd_dir_modified, spatial_dimension, p, para)
    output = table('Size',[1 18], ...
        'VariableTypes',["cell", "cell", "cell", "cell", "double", "double", ...
        "cell", "double", "cell", "double", "double", "double", ....
        "cell", "cell", "double", "cell", "cell", "double"], ...
        'VariableNames',["map", "fields_map", "fields", "map_stability", "information_content", "information_rate", ...
        "autocorrelogram", "grid_score", "grid_stat", "center_field", "score_radius", "border_score", ...
        "tc", "tcStat", "curMVL", "angle_stability", "speed_tuning_curve", "speed_score"]);
    if do_cell_analysis.place_cell || do_cell_analysis.grid_cell || do_cell_analysis.border_cell
        spike_pos_temp = spike_pos(:,1:3);
        if spatial_dimension(1) == 2
            spike_pos_temp(isnan(spike_pos_temp(:,2)),:) = [];
            if ~isempty(spike_pos_temp)
                % calculate rate map;
                output.map{1} = analyses.map(pos_modified, spike_pos_temp(:,1), para.map{:});
                output.map{1}.spike_pos = spike_pos_temp;
                % 'datatime', p.datatime, 'binWidth', p.binWidth, 'smooth', p.smooth, ...
                % 'minTime', p.binMinTime, 'blanks', p.blank, 'maxGap', p.maxGap, 'limits', p.limits
                
                % field;
                [output.fields_map{1}, output.fields{1}] = analyses.placefield(output.map{1}, 'pos',pos_modified, para.field{:});
                % 'binWidth', p.binWidth, 'pos', pos, 'minPeak', 0.5, 'threshold', 0.6
                
                % stability;
                output.map_stability{1} = cellStability(spike_pos_temp, splitDataArray_modified, para.map);
            end
        elseif spatial_dimension(1) == 1
            spike_pos_temp(isnan(spike_pos_temp(:,spatial_dimension(2)+1)),:) = [];
            if ~isempty(spike_pos_temp)
                % calculate rate map;
                output.map{1} = analyses.map(pos_modified, spike_pos_temp(:,1), para.map{:});
                output.map{1}.spike_pos = spike_pos_temp;
                
                % field;
                [output.fields_map{1}, output.fields{1}] = analyses.placefield1D(output.map{1}, 'pos',pos_modified, para.field{:});
                % 'binWidth', p.binWidth, 'pos', pos  %% attention the difference of input in placefield!
                
                % stability;
                output.map_stability{1} = cellStability(spike_pos_temp(:,[1,spatial_dimension(2)+1]), splitDataArray_modified, para.map);
            end
        end
    end
    
    if do_cell_analysis.place_cell
        if ~isempty(spike_pos_temp)
            % information;
            [information, ~, ~] = analyses.mapStatsPDF(output.map{1});
            output.information_content = information.content;
            output.information_rate = information.rate;
        end
        
    end
    
    if do_cell_analysis.grid_cell
        if ~isempty(spike_pos_temp)
            % 2D autocorrelation (autocorrelogram) of a firing map;
            output.autocorrelogram{1} = analyses.autocorrelation(output.map{1}.z);
            
            % gridness score;
            [output.grid_score, output.grid_stat{1}, output.center_field, ~, output.score_radius] = ...
                analyses.gridnessScore(output.autocorrelogram{1}, para.grid{:});
        end
    end
    
    if do_cell_analysis.border_cell
        if ~isempty(spike_pos_temp)
            output.border_score = analyses.borderScore(output.map{1}.z, output.fields_map{1}, output.fields{1}, para.border{:});
        end
    end
    
    if do_cell_analysis.head_direction_cell
        spike_pos_hd = spike_pos(:, [1,4]);
        spike_pos_hd(isnan(spike_pos_hd(:, 2)),:) = [];
        if isempty(spike_pos_hd) || isempty(hd_dir_modified)
        else
            % turning curve;
            output.tc{1} = analyses.turningCurve(spike_pos_hd(:, 2), hd_dir_modified, p.hd2.sampleTime, para.hd{:});
            output.tcStat{1} = analyses.tcStatistics(output.tc{1}, p.hd.binWidth, p.hd2.percentile);
            output.tcStat{1}.spkDir = spike_pos_hd;
            output.tcStat{1}.p.hd = p.hd;
            output.tcStat{1}.p.hd2 = p.hd2;
            
            % mean vector length;
            output.curMVL = output.tcStat{1}.r;
            
            % stability;
            output.angle_stability{1} = cellStability_hd(spike_pos_hd, splitDataArray_hd_modified, p, para);
        end
    end
    
    if do_cell_analysis.speed_cell
        speed_firing_rate = p.speed2.speed_firing_rate(:,cell_i);
        output.speed_score = corr(speed_firing_rate, p.speed2.speed_input, 'rows','pairwise');
        % analyses.speedScore(speed_modified, p.speed2.firing_rate(:,cell_i), p.speed2.span, para.speed{:});
        speed_range = p.speed2.speed_range_bin;
        speed_center = (speed_range(1:end-1) + speed_range(2:end)) / 2;
        speed_smooth_bin = floor((p.speed2.speed_input - speed_range(1)) / p.speed2.speed_bin);
        speed_tuning_curve = cell2mat(arrayfun(@(x) nanmean(speed_firing_rate(speed_smooth_bin == x)), ...
            (speed_range(1:end-1)' - speed_range(1)) / p.speed2.speed_bin, 'UniformOutput',0));
        output.speed_tuning_curve{1} = [speed_center' speed_tuning_curve];
    end
end

function output_shuffle = cellShuffle(output, do_cell_analysis, spike_pos_shuffle, ...
        splitDataArray_modified, splitDataArray_hd_modified, pos_modified_full, hd_dir_modified, ...
        spatial_dimension, fr_shuffle, speed_input, p, para)
    output_shuffle = table('Size',[1 10], ...
        'VariableTypes',["double", "double", "double", "double", ...
        "double", "double", "double", "double", "double", "double"], ...
        'VariableNames',["map_stability_half", "map_stability_time", ...
        "information_content", "information_rate", "grid_score", "border_score", ...
        "curMVL", "angle_stability_half", "angle_stability_time", "speed_score"]);
    % "spike_pos", "map", "fields_map", "fields",
    
    if do_cell_analysis.place_cell || do_cell_analysis.grid_cell || do_cell_analysis.border_cell
        if spatial_dimension(1) == 2
            pos_modified = pos_modified_full;
            spike_pos_temp = spike_pos_shuffle(:,1:3);
        elseif spatial_dimension(1) == 1
            pos_modified = pos_modified_full(:,[1,spatial_dimension(2)+1]);
            spike_pos_temp = spike_pos_shuffle(:,[1,spatial_dimension(2)+1]);
        end
        
        spike_pos_temp(isnan(spike_pos_temp(:,2)),:) = [];
        if ~isempty(spike_pos_temp)
            % calculate rate map;
            map = analyses.map(pos_modified, spike_pos_temp(:,1), para.map{:});
            
            if spatial_dimension(1) == 2
                % field;
                [fields_map, fields] = analyses.placefield(map, 'pos', pos_modified, para.field{:});
                % 'binWidth', p.binWidth, 'pos', pos, 'minPeak', 0.5, 'threshold', 0.6
            end
            
            % stability;
            map_stability_shuffle = cellStability(spike_pos_temp, splitDataArray_modified, para.map);
            output_shuffle.map_stability_half = map_stability_shuffle.half_correlation;
            output_shuffle.map_stability_time = map_stability_shuffle.time_correlation;
        end
    end
    
    if do_cell_analysis.place_cell
        if ~isempty(spike_pos_temp)
            % information;
            [information, ~, ~] = analyses.mapStatsPDF(map);
            output_shuffle.information_content = information.content;
            output_shuffle.information_rate = information.rate;
        end
    end
    
    if do_cell_analysis.grid_cell
        if ~isempty(spike_pos_temp)
            % 2D autocorrelation (autocorrelogram) of a firing map;
            autocorrelogram_shuffle = analyses.autocorrelation(map.z);
            
            % gridness score;
            output_shuffle.grid_score = analyses.gridnessScoreShuffled(autocorrelogram_shuffle, ...
                output.center_field, output.score_radius, p.grid2.radii);
        end
    end
    
    if do_cell_analysis.border_cell
        if ~isempty(spike_pos_temp)
            output_shuffle.border_score = ...
                analyses.borderScore(map.z, fields_map, fields, para.border{:});
        end
    end
    
    if do_cell_analysis.head_direction_cell
        spike_pos_hd = spike_pos_shuffle(:,[1,4]);
        spike_pos_hd(isnan(spike_pos_hd(:,2)),:) = [];
        if isempty(spike_pos_hd) || isempty(hd_dir_modified)
        else
            % turning curve;
            tc = analyses.turningCurve(spike_pos_hd(:,2), hd_dir_modified, p.hd2.sampleTime, para.hd{:});
            tcStat = analyses.tcStatistics(tc, p.hd.binWidth, p.hd2.percentile);
            
            % mean vector length;
            output_shuffle.curMVL = tcStat.r;
            
            % stability;
            angle_stability_shuffle = cellStability_hd(spike_pos_hd, splitDataArray_hd_modified, p, para);
            output_shuffle.angle_stability_half = angle_stability_shuffle.half_correlation;
            output_shuffle.angle_stability_time = angle_stability_shuffle.time_correlation;
        end
    end
    
    if do_cell_analysis.speed_cell
        output_shuffle.speed_score = corr(fr_shuffle, speed_input, 'rows','pairwise');
        % speed_score = analyses.speedScore(speed_modified, p.speed2.firing_rate(:,cell_i), p.speed2.span, para.speed{:});
        % output_shuffle.speed_score = speed_score(2);
    end
end

function output_shuffle = shuffleProcess(do_shuffle, do_cell_analysis, cell_i, k, calcium_time, ...
        shuffle_num, threshold_prc, output_temp, pos_modified_full, hd_dir_modified, spatial_dimension, ...
        splitDataArray_modified, splitDataArray_hd_modified, p, para)
    output_shuffle = table('Size',[shuffle_num + 3, 10], ...
        'VariableTypes',["double", "double", "double", "double", ...
        "double", "double", "double", "double", "double", "double"], ...
        'VariableNames',["map_stability_half", "map_stability_time", ...
        "information_content", "information_rate", "grid_score", "border_score", ...
        "curMVL", "angle_stability_half", "angle_stability_time", "speed_score"]);
    if ~do_shuffle, output_shuffle{end, output_shuffle{end,:} == 0} = nan; return; end
    
    time_30s = find(calcium_time >= 30 + calcium_time(1), 1); % 30 s;
    %dis_shuffle = cell(shuffle_num, 1);
    time_shuffle_i = time_30s + randi((length(calcium_time) - time_30s*2), shuffle_num,1);
    
    for m = 1:shuffle_num
        % shuffle calcium time;
        k_shuffle = circshift(k, time_shuffle_i(m));
        if do_cell_analysis.speed_cell, fr_shuffle = circshift(p.speed2.speed_firing_rate(:,cell_i), time_shuffle_i(m));
        else, fr_shuffle = []; p.speed2.speed_input = []; end
        spike_pos_shuffle = calcium_time(k_shuffle); % mod(k + time_shuffle_i(m), length(calcium_time))
        spike_pos_shuffle = spikePos(spike_pos_shuffle, pos_modified_full(:,1), pos_modified_full, hd_dir_modified);
        if isempty(spike_pos_shuffle), continue; end
        
        output_shuffle(m, :) = cellShuffle(output_temp, do_cell_analysis, spike_pos_shuffle, ...
            splitDataArray_modified, splitDataArray_hd_modified, pos_modified_full, hd_dir_modified, ...
            spatial_dimension, fr_shuffle, p.speed2.speed_input, p, para);
    end
    
    % 95% and 99% threshold;
    output_shuffle = calcThreshold(output_shuffle, shuffle_num, threshold_prc);
end

function [cell_all, output_shuffle, cell_sig] = cellDetectFcn(do_cell_analysis, cell_i, output_temp, output_shuffle_temp)
    cell_all = table('Size',[1, 5], ...
        'VariableTypes',["double", "double", "double", "double", "double"], ...
        'VariableNames',["place_cell", "grid_cell", "border_cell", "head_direction_cell", "speed_cell"]);
    cell_sig = struct;
    output_shuffle = table('Size',[1, 10], ...
        'VariableTypes',["double", "double", "double", "double", ...
        "double", "double", "double", "double", "double", "double"], ...
        'VariableNames',["map_stability_half_shuffle", "map_stability_time_shuffle", ...
        "information_content_shuffle", "information_rate_shuffle", "grid_score_shuffle", "border_score_shuffle", ...
        "curMVL_shuffle", "angle_stability_half_shuffle", "angle_stability_time_shuffle", "speed_score_shuffle"]);
    
    if do_cell_analysis.place_cell || do_cell_analysis.grid_cell || do_cell_analysis.border_cell
        output_shuffle.map_stability_half_shuffle = output_shuffle_temp.map_stability_half';
        output_shuffle.map_stability_time_shuffle = output_shuffle_temp.map_stability_time';
    end
    
    if do_cell_analysis.place_cell
        % information_all{i,:} = [output_temp.information_content, output_temp.information_rate];
        output_shuffle.information_content_shuffle = output_shuffle_temp.information_content';
        output_shuffle.information_rate_shuffle = output_shuffle_temp.information_rate';
        if output_temp.information_content > output_shuffle_temp.information_content(end) % && map_stability.half_correlation > 0.3
            % || map_stability.time_correlation > 0.3)
            cell_all.place_cell = cell_i;
            cell_sig.place_cell = 1;
        else
            cell_sig.place_cell = 0;
        end
    end
    
    if do_cell_analysis.grid_cell
        output_shuffle.grid_score_shuffle = output_shuffle_temp.grid_score';
        if output_temp.grid_score > output_shuffle_temp.grid_score(end)
            cell_all.grid_cell = cell_i;
            cell_sig.grid_cell = 1;
        else
            cell_sig.grid_cell = 0;
        end
    end
    
    if do_cell_analysis.border_cell
        output_shuffle.border_score_shuffle = output_shuffle_temp.border_score';
        if output_temp.border_score > output_shuffle_temp.border_score(end)
            cell_all.border_cell = cell_i;
            cell_sig.border_cell = 1;
        else
            cell_sig.border_cell = 0;
        end
    end
    
    if do_cell_analysis.head_direction_cell
        output_shuffle.curMVL_shuffle = output_shuffle_temp.curMVL';
        output_shuffle.angle_stability_half_shuffle = output_shuffle_temp.angle_stability_half';
        output_shuffle.angle_stability_time_shuffle = output_shuffle_temp.angle_stability_time';
        
        if output_temp.curMVL > output_shuffle_temp.curMVL(end)
            cell_all.head_direction_cell = cell_i;
            cell_sig.head_direction_cell = 1;
        else
            cell_sig.head_direction_cell = 0;
        end
    end
    
    if do_cell_analysis.speed_cell
        output_shuffle.speed_score_shuffle = output_shuffle_temp.speed_score';
        if output_temp.speed_score > output_shuffle_temp.speed_score(end)
            cell_all.speed_cell = cell_i;
            cell_sig.speed_cell = 1;
        else
            cell_sig.speed_cell = 0;
        end
    end
end

function output_shuffle = calcThreshold(output_shuffle, shuffle_num, threshold_prc)
    for j = 1:size(output_shuffle,2)
        output_type = output_shuffle.Properties.VariableNames{j};
        output_shuffle{output_shuffle{:,j} == 0, j} = nan;
        output_shuffle{shuffle_num + 1, j} = prctile(output_shuffle{1:shuffle_num, j}, 95);
        output_shuffle{shuffle_num + 2, j} = prctile(output_shuffle{1:shuffle_num, j}, 99);
        output_shuffle{shuffle_num + 3, j} = prctile(output_shuffle{1:shuffle_num, j}, threshold_prc.(output_type));
    end
end

%-------------------------------------------------------------------------%
%                                 plots                                   %
%-------------------------------------------------------------------------%
% cell;
function singleCellPlot(SFP, cell_i, save_folder, fig_fmt, color_cell)
    markPartSFP_color_v3(SFP, SFP(cell_i,:,:), color_cell);
    for fmt_i = 1:length(fig_fmt)
        saveas(gcf, strcat(save_folder, '/cell/cell', num2str(cell_i), 'inAll', fig_fmt{fmt_i}));
    end
end

% map;
function ratemapPlot(cell_i, spatial_dimension, map, save_folder, fig_fmt)
    figure;
    if spatial_dimension(1) == 1
        if spatial_dimension(2) == 1, plot.colorMap(map.z);
        elseif spatial_dimension(2) == 2,  plot.colorMap(map.z'); end
    else
        plot.colorMap(map.z);
    end
    axis equal;
    axis off;
    title(sprintf('%s%3.2f','Rate Map. Peak = ',nanmax(nanmax(map.z))));
    for fmt_i = 1:length(fig_fmt)
        saveas(gcf,strcat(save_folder, '/ratemap/cell', num2str(cell_i), 'ratemap', fig_fmt{fmt_i}));
    end
    
    if spatial_dimension(1) == 1
        figure; plot(map.z, 'Color',[0.3 0.3 0.3]);
        xlabel('Distance (cm)');
        ylabel('Firing rate (Hz)');
        box off;
        for fmt_i = 1:length(fig_fmt)
            saveas(gcf,strcat(save_folder, '/ratemap/cell', num2str(cell_i), 'ratetuning', fig_fmt{fmt_i}));
        end
    end
end

function cellpathPlot(cell_i, pos, behavLimit, spike_pos, save_folder, fig_fmt)
    figure;
    plot(pos(:,2),pos(:,3), 'color',[0.5,0.5,0.5],'LineWidth',1);
    hold on;
    scatter(spike_pos(:,2),spike_pos(:,3), 20, 'r', 'filled');
    hold off;
    axis equal;
    xlimit = get(gca,'xlim');
    ylimit = get(gca,'ylim');
    if xlimit(2)-xlimit(1) < behavLimit(2) - behavLimit(1), xlim(behavLimit(1:2)); end
    if ylimit(2)-ylimit(1) < behavLimit(4) - behavLimit(3), ylim(behavLimit(3:4)); end
    % for fmt_i = 1:length(fig_fmt)
    %     saveas(gcf,strcat(save_folder, '/cellpath/cell', num2str(cell_i), 'path', fig_fmt{fmt_i})); % .fig
    % end
    
    axis off;
    for fmt_i = 1:length(fig_fmt)
        saveas(gcf,strcat(save_folder, '/cellpath/cell', num2str(cell_i), 'path-noAxis', fig_fmt{fmt_i}));
    end
end

% grid cell;
function autocorrelogramPlot(cell_i, autocorrelogram, grid_score, save_folder, save_subfolder, fig_fmt)
    figure;
    plot.colorMap(autocorrelogram, 'cutoffs',[-1, 1]);
    axis equal;
    axis off;
    title(sprintf('%s%3.2f', 'Autocorrelogram. Grid = ', grid_score));
    for fmt_i = 1:length(fig_fmt)
        saveas(gcf,strcat(save_folder, '/', save_subfolder, ...
            '/autocorrelogram/cell', num2str(cell_i), 'autocorrelogram', fig_fmt{fmt_i}));
    end
end

% head direction cell;
function turningcurvePlot(cell_i, tc, tcStat, cell_sig, save_folder, save_subfolder, fig_fmt, p)
    curveToDisplay = tc(:, 2) / max(tc(:, 2));
    
    % figure;
    figure;
    hold on;
    plot.circularTurning(curveToDisplay, '-k');
    plot.circularTurning(p.trajectory_hd, p.trajectoryPlot{:});
    hold off;
    if cell_sig.head_direction_cell
        title(sprintf('Peak: %3.2f, MVL: %3.2f*', tcStat.peakRate, tcStat.r));
    else
        title(sprintf('Peak: %3.2f, MVL: %3.2f', tcStat.peakRate, tcStat.r));
    end
    axis equal;
    
    for fmt_i = 1:length(fig_fmt)
        saveas(gcf,strcat(save_folder, '/', save_subfolder, ...
            '/turningcurve/cell', num2str(cell_i), 'turningcurve', fig_fmt{fmt_i}));
    end
end

function cellpathPlot_hd(cell_i, pos, behavLimit, spike_pos, save_folder, save_subfolder, fig_fmt)
    figure;
    plot(pos(:,2),pos(:,3), 'color',[0.5,0.5,0.5], 'LineWidth',1);
    hold on;
    scatter(spike_pos(:,2),spike_pos(:,3), 20, spike_pos(:,4), 'filled');
    colormap hsv; caxis([0 360]);
    hold off;
    axis equal;
    xlimit = get(gca,'xlim');
    ylimit = get(gca,'ylim');
    if xlimit(2)-xlimit(1) < behavLimit(2) - behavLimit(1), xlim(behavLimit(1:2)); end
    if ylimit(2)-ylimit(1) < behavLimit(4) - behavLimit(3), ylim(behavLimit(3:4)); end
    axis off;
    
    for fmt_i = 1:length(fig_fmt)
        saveas(gcf,strcat(save_folder, '/', save_subfolder, ...
            '/cellpath/cell', num2str(cell_i), 'path-noAxis', fig_fmt{fmt_i}));
    end
end

% merge;
function plotMerge(do_cell_analysis, spatial_dimension, cell_i, output, cell_sig, pos, p, save_folder, fig_fmt)
    figure('Position', [600 250 700 500]);
    
    subplot('Position', [0.05 0.52 0.4 0.4]);
    if spatial_dimension(1) == 1
        if spatial_dimension(2) == 1, plot.colorMap(output.map{1}.z, 'bar','off');
        elseif spatial_dimension(2) == 2,  plot.colorMap(output.map{1}.z', 'bar','off'); end
    else, plot.colorMap(output.map{1}.z, 'bar','off');
    end
    axis equal;
    axis off;
    if do_cell_analysis.place_cell
        title_text1 = sprintf(', Info: %3.2f', output.information_content);
        if cell_sig.place_cell, title_text1 = [title_text1, '*']; end
    else
        title_text1 = [];
    end
    if do_cell_analysis.border_cell
        title_text2 = sprintf(', Border: %3.2f', output.border_score);
        if cell_sig.border_cell, title_text2 = [title_text2, '*']; end
    else
        title_text2 = [];
    end
    title(sprintf('Peak: %3.2f%s%s', nanmax(nanmax(output.map{1}.z)), title_text1, title_text2));
    
    if do_cell_analysis.grid_cell
        subplot('Position', [0.55 0.52 0.4 0.4]);
        plot.colorMap(output.autocorrelogram{1}, 'bar','off', 'cutoffs',[-1, 1]);
        axis equal;
        axis off;
        if cell_sig.grid_cell, title(sprintf('Grid: %3.2f*', output.grid_score));
        else, title(sprintf('Grid: %3.2f', output.grid_score));
        end
    end
    
    subplot('Position', [0.05 0.05 0.4 0.4]);
    plot(pos(:,2),pos(:,3), 'color',[0.5,0.5,0.5], 'LineWidth',1);
    hold on;
    scatter(output.map{1}.spike_pos(:,2),output.map{1}.spike_pos(:,3), 10, 'r', 'filled');
    % 'MarkerFaceColor',[1 0 0], 'MarkerEdgeColor',[1 0 0]
    hold off;
    axis equal;
    xlimit = get(gca,'xlim');
    ylimit = get(gca,'ylim');
    if xlimit(2)-xlimit(1) < p.behav_limit(2) - p.behav_limit(1), xlim(p.behav_limit(1:2)); end
    if ylimit(2)-ylimit(1) < p.behav_limit(4) - p.behav_limit(3), ylim(p.behav_limit(3:4)); end
    axis off;
    
    if do_cell_analysis.head_direction_cell
        curveToDisplay = output.tc{1}(:, 2) / max(output.tc{1}(:, 2));
        subplot('Position', [0.55 0.05 0.4 0.4]);
        hold on;
        plot.circularTurning(curveToDisplay, '-k');
        plot.circularTurning(p.hd2.trajectory_hd, p.hd2.trajectoryPlot{:});
        hold off;
        if cell_sig.head_direction_cell
            title(sprintf('Peak: %3.2f, MVL: %3.2f*', output.tcStat{1}.peakRate, output.tcStat{1}.r));
        else
            title(sprintf('Peak: %3.2f, MVL: %3.2f', output.tcStat{1}.peakRate, output.tcStat{1}.r));
        end
        axis equal;
    end
    
    for fmt_i = 1:length(fig_fmt)
        saveas(gcf,strcat(save_folder, '/merge/cell', num2str(cell_i), fig_fmt{fmt_i}));
    end
end

% speed tuning;
function speedTuningPlot(cell_i, speed_tuning_curve, ...
        speed_score, cell_sig, save_folder, save_subfolder, fig_fmt)
    figure; plot(speed_tuning_curve(:,1), speed_tuning_curve(:,2)); box off;
    xlabel('Speed (cm/s)'); ylabel('Firing rate (Hz)');
    if cell_sig.speed_cell, title(['Score: ', num2str(speed_score), '*']);
    else, title(['Score: ', num2str(speed_score)]); end
    
    for fmt_i = 1:length(fig_fmt)
        saveas(gcf,strcat(save_folder, '/', save_subfolder, ...
            '/speedtuning/cell', num2str(cell_i), 'speedtuning', fig_fmt{fmt_i}));
    end
end

function speedTuningPlot_all(speed_firing_rate, speed_input, speed_range, speed_bin, save_folder, save_subfolder, fig_fmt)
    % speed_bin = 2;
    % speed_range = 0:speed_bin:max(speed_input);
    % speed_range = [speed_range, speed_range(end)+speed_bin];
    speed_center = (speed_range(1:end-1) + speed_range(2:end)) / 2;
    speed_smooth_bin = floor((speed_input - speed_range(1)) / speed_bin);
    speed_tuning_curve = cell2mat(arrayfun(@(n) arrayfun(@(x) nanmean(speed_firing_rate(speed_smooth_bin == x, n)), ...
        (speed_range(1:end-1)' - speed_range(1)) / speed_bin), 1:size(speed_firing_rate,2), 'UniformOutput',0));
    
    figure; plot(speed_center, speed_tuning_curve'); box off;
    xlabel('Speed (cm/s)'); ylabel('Firing rate (Hz)');
    for fmt_i = 1:length(fig_fmt)
        saveas(gcf,strcat(save_folder, '/', save_subfolder, ...
            '/speedtuning/cell_speedtuning_all', fig_fmt{fmt_i}));
    end
    
    figure; plot(speed_center, nanmean(speed_tuning_curve,2)); box off;
    xlabel('Speed (cm/s)'); ylabel('Firing rate (Hz)');
    for fmt_i = 1:length(fig_fmt)
        saveas(gcf,strcat(save_folder, '/', save_subfolder, ...
            '/speedtuning/cell_speedtuning_mean', fig_fmt{fmt_i}));
    end
    
end

function cellPropPlot(draw_fig, do_cell_analysis, spatial_dimension, cell_i, output_temp, SFP, pos, spike_pos, ...
        cell_sig, cell_name, color_cell, p, save_folder, fig_fmt)
    if draw_fig.cell, singleCellPlot(SFP, cell_i, save_folder, fig_fmt, color_cell.cell); end
    
    if do_cell_analysis.place_cell || do_cell_analysis.grid_cell || do_cell_analysis.border_cell
        if draw_fig.ratemap, ratemapPlot(cell_i, spatial_dimension, output_temp.map{1}, save_folder, fig_fmt); end
        if draw_fig.cellpath, cellpathPlot(cell_i, pos, p.behav_limit, spike_pos, save_folder, fig_fmt); end
        if draw_fig.merge, plotMerge(do_cell_analysis, spatial_dimension, cell_i, output_temp, cell_sig, ...
                pos, p, save_folder, fig_fmt); end
    end
    
    if do_cell_analysis.place_cell
    end
    
    if do_cell_analysis.grid_cell
        if draw_fig.autocorrelogram, autocorrelogramPlot(cell_i, output_temp.autocorrelogram{1}, ...
                output_temp.grid_score, save_folder, cell_name.grid_cell, fig_fmt); end
    end
    
    if do_cell_analysis.border_cell
    end
    
    if do_cell_analysis.head_direction_cell
        if draw_fig.turningcurve, turningcurvePlot(cell_i, output_temp.tc{1}, output_temp.tcStat{1}, ...
                cell_sig, save_folder, cell_name.head_direction_cell, fig_fmt, p.hd2);end
        if draw_fig.cellpath_hd, cellpathPlot_hd(cell_i, pos, p.behav_limit, spike_pos, ...
                save_folder, cell_name.head_direction_cell, fig_fmt); end
    end
    
    if do_cell_analysis.speed_cell
        if draw_fig.speedtuning
            speedTuningPlot(cell_i, output_temp.speed_tuning_curve{1}, ...
                output_temp.speed_score, cell_sig, save_folder, cell_name.speed_cell, fig_fmt);
        end
    end
    
    close all;
end

function cellsPlot(cell_all, SFP_all, save_folder, save_subfolder, fig_fmt, color_cell)
    if isempty(cell_all), return; end
    SFP = SFP_all(cell_all,:,:);
    save(strcat(save_folder, '/', save_subfolder, '/SFP.mat'), 'SFP');
    
    markPartSFP_color_v3(SFP_all, SFP, color_cell);
    for fmt_i = 1:length(fig_fmt)
        saveas(gcf,strcat(save_folder, '/', save_subfolder, '/cellsinAll', fig_fmt{fmt_i}));
    end
    close();
end

%-------------------------------------------------------------------------%
%                            output process                               %
%-------------------------------------------------------------------------%
function cell_all = saveCell(do_shuffle, cell_all, cell_name, save_folder)
    if ~do_shuffle, cell_all = []; return; end
    cell_all(isnan(cell_all)) = [];
    cell_all(cell_all == 0) = [];
    eval([cell_name ' = cell_all;']);
    save([save_folder, '/', cell_name, '/', cell_name, '.mat'], cell_name);
end

function saveInfo(output, info_names, save_folder, save_subfolder) %#ok<INUSL>
    for j = 1:length(info_names)
        info_name = info_names{j};
        eval([info_name, ' = output.', info_name, ';']);
        save([save_folder, '/', save_subfolder, '/', info_name, '.mat'], info_name);
    end
end

function saveShuffle(do_shuffle, output_shuffle, shuffle_names, save_folder, save_subfolder) %#ok<INUSL>
    if ~do_shuffle, return; end
    for j = 1:length(shuffle_names)
        info_name = shuffle_names{j};
        eval([info_name, ' = output_shuffle{j};']);
        save([save_folder, '/', save_subfolder, '/', info_name, '.mat'], info_name);
    end
end

function saveCellResults(do_cell_analysis, do_shuffle, output, output_shuffle, cell_all, cell_name, ...
        save_folder, draw_fig, SFP, fig_fmt, color_cell, p)
    if do_cell_analysis.place_cell || do_cell_analysis.grid_cell || do_cell_analysis.border_cell
        saveInfo(output, {'map', 'fields_map', 'fields', 'map_stability'}, save_folder, 'map');
        saveShuffle(do_shuffle, ...
            {output_shuffle.map_stability_half_shuffle', output_shuffle.map_stability_time_shuffle'}, ...
            {'map_stability_half_shuffle', 'map_stability_time_shuffle'}, save_folder, 'map');
    end
    
    if do_cell_analysis.place_cell
        place_cell = saveCell(do_shuffle, cell_all.place_cell, cell_name.place_cell, save_folder);
        
        saveInfo(output, {'information_content', 'information_rate'}, save_folder, cell_name.place_cell);
        saveShuffle(do_shuffle, ...
            {output_shuffle.information_content_shuffle', output_shuffle.information_rate_shuffle'}, ...
            {'information_content_shuffle', 'information_rate_shuffle'}, save_folder, cell_name.place_cell);
        
        if draw_fig.cell, cellsPlot(place_cell, SFP, save_folder, cell_name.place_cell, fig_fmt, color_cell.place_cell); end
    end
    
    if do_cell_analysis.grid_cell
        grid_cell = saveCell(do_shuffle, cell_all.grid_cell, cell_name.grid_cell, save_folder);
        
        saveInfo(output, {'autocorrelogram', 'grid_score', 'grid_stat', 'center_field', 'score_radius'}, ...
            save_folder, cell_name.grid_cell);
        saveShuffle(do_shuffle, {output_shuffle.grid_score_shuffle'}, ...
            {'grid_score_shuffle'}, save_folder, cell_name.grid_cell);
        
        if draw_fig.cell, cellsPlot(grid_cell, SFP, save_folder, cell_name.grid_cell, fig_fmt, color_cell.grid_cell); end
    end
    
    if do_cell_analysis.border_cell
        border_cell = saveCell(do_shuffle, cell_all.border_cell, cell_name.border_cell, save_folder);
        
        saveInfo(output, {'border_score'}, save_folder, cell_name.border_cell);
        saveShuffle(do_shuffle, {output_shuffle.border_score_shuffle'}, ...
            {'border_score_shuffle'}, save_folder, cell_name.border_cell);
        
        if draw_fig.cell, cellsPlot(border_cell, SFP, save_folder, cell_name.border_cell, ...
                fig_fmt, color_cell.border_cell); end
    end
    
    if do_cell_analysis.head_direction_cell
        head_direction_cell = saveCell(do_shuffle, cell_all.head_direction_cell, cell_name.head_direction_cell, save_folder);
        
        saveInfo(output, {'curMVL', 'tc', 'tcStat', 'angle_stability'}, save_folder, cell_name.head_direction_cell);
        saveShuffle(do_shuffle, {output_shuffle.curMVL_shuffle', ...
            output_shuffle.angle_stability_half_shuffle', output_shuffle.angle_stability_time_shuffle'}, ...
            {'curMVL_shuffle', 'angle_stability_half_shuffle', 'angle_stability_time_shuffle'}, ...
            save_folder, cell_name.head_direction_cell);
        trajectory_map = p.hd2.trajectory_map;
        save([save_folder, '/', cell_name.head_direction_cell, '/trajectory_map.mat'], 'trajectory_map');
        
        if draw_fig.cell, cellsPlot(head_direction_cell, SFP, save_folder, cell_name.head_direction_cell, ...
                fig_fmt, color_cell.head_direction_cell); end
    end
    
    if do_cell_analysis.speed_cell
        speed_cell = saveCell(do_shuffle, cell_all.speed_cell, cell_name.speed_cell, save_folder);
        
        saveInfo(output, {'speed_score'}, save_folder, cell_name.speed_cell);
        saveShuffle(do_shuffle, {output_shuffle.speed_score_shuffle'}, ...
            {'speed_score_shuffle'}, save_folder, cell_name.speed_cell);
        speed_firing_rate = p.speed2.speed_firing_rate;
        speed_input = p.speed2.speed_input;
        save([save_folder, '/', cell_name.speed_cell, '/speed_firing_rate.mat'], 'speed_firing_rate');
        save([save_folder, '/', cell_name.speed_cell, '/speed_input.mat'], 'speed_input');
        
        if draw_fig.cell, cellsPlot(speed_cell, SFP, save_folder, cell_name.speed_cell, ...
                fig_fmt, color_cell.speed_cell); end
        if draw_fig.speedtuning, speedTuningPlot_all(speed_firing_rate, speed_input, ...
                p.speed2.speed_range_bin, p.speed2.speed_bin, save_folder, cell_name.speed_cell, fig_fmt); end
    end
end
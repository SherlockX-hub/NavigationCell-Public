%% cellDecoding.m
% Decode positions with all cells.
% Use GNBDecoder written by Kai Gao.

% v3: May, 2023.
% Corss-validation version.

% Created by Xiang Zhang, 2022.
clear;

%% code path;
% addpath('D:\Calcuim_Imaing_Code\GK_Inscopix_Auto_PCA-ICA\CAN');
% addpath(genpath('F:\LJY\retro-MEC-HPC\BNT'));

%% parameters;
cell_type1 = 'self_cell';
cell_type2 = 'place_cell';
% training_set = [0 0.8];%[0 0.5; 0.65 1];
% testing_set = [0.8 1];
bin_size = 0.2;
PC_num = 1;

p_map.binWidth = 2;
p_map.smooth = 2;
p_map.minTime = 0;
p_map.maxGap = 0.3;

[fit_loss_PC, predict_loss_PC] = deal(cell(1, length(cell_type1)));
[fit_loss_PC_shuffle, predict_loss_PC_shuffle] = deal(cell(1, length(cell_type1)));

%% main function;
session_num = 1;
session_path = 'G:\ZX\Data_temp\49-20230108\49-20230108-1';
sInd = strfind(session_path, '\');
session_ID = session_path(sInd(end-1):end);

load([session_path, '\NeuronActivity.mat'], 'NeuronActivity');
load([session_path, '\behav.mat'], 'behav');

if NeuronActivity.timestamps(1) - behav.timestamps(1) == 0
    calcium_time = NeuronActivity.time;
    behav_time = behav.time;
elseif NeuronActivity.timestamps(1) > behav.timestamps(1)
    calcium_time = seconds(NeuronActivity.timestamps - behav.timestamps(1));
    behav_time = behav.time;
elseif NeuronActivity.timestamps(1) < behav.timestamps(1)
    calcium_time = NeuronActivity.time;
    behav_time = seconds(behav.timestamps - NeuronActivity.timestamps(1));
end
calcium_event = NeuronActivity.Event_filtered_exp2;

p_map.limits = [0 behav.trackLength(1) 0 behav.trackLength(2)];

[~, pos] = loadInfo(cell_type1, cell_type2, session_path); % PC
cell_sig = 1:size(calcium_event,2);% all cells;
% if length(cell_sig) < PC_num, continue; end
% if isempty(pos), continue; end

% spike time stamp;
spike_timestamps_PC = findPCSpikeTime(cell_sig, calcium_time, calcium_event);
% if isempty(spike_timestamps_PC), continue; end

% decode;
[fit_loss_PC{1}, predict_loss_PC{1}] = ...
    decodeIntegration(spike_timestamps_PC, pos, p_map, ...
    bin_size, 0, [num2str(session_num), '-', session_ID], [cell_type1, '-', cell_type2]);

% shuffle;
[fit_loss_PC_shuffle{1}, predict_loss_PC_shuffle{1}] = ...
    decodeShuffle(spike_timestamps_PC, pos, p_map, bin_size);

% done;


%% save files;
save('pcover1_cv.mat', ...
    'fit_loss_PC','fit_loss_PC_shuffle','predict_loss_PC','predict_loss_PC_shuffle');

%% functions;
function [cell_sig, pos] = loadInfo(cell_type1, cell_type2, session_temp)
    try
        load([session_temp, filesep, cell_type1, filesep, cell_type2, filesep, cell_type2, '.mat'], cell_type2);
        load([session_temp, filesep, cell_type1, filesep, 'pos.mat'], 'pos');
        eval(['cell_sig = ', cell_type2, ';']);
    catch
        cell_sig = [];
        pos = nan(0,1);
    end
end

function spike_timestamps = findPCSpikeTime(cell_sig, calcium_time, calcium_event)
    if isempty(cell_sig), spike_timestamps = {}; return; end
    
    spike_timestamps = cell(1, length(cell_sig));
    for PC_i = 1:length(cell_sig)
        i = cell_sig(PC_i);
        k = find(calcium_event(:,i) > 0); % 3 * std(calcium_event(:,i))
        if ~isempty(k) && length(k) > 10
            spike_timestamps{PC_i} = calcium_time(k); %pos(knnsearch(pos(:, 1), ms.time(k)), 1);
        end
    end
    spike_timestamps(cellfun(@isempty, spike_timestamps(:))) = [];
end

function [spike_set_PC, pos_set] = divideData(spike_timestamps, pos, time_set)
    spike_set_PC = cell(length(spike_timestamps), 1);
    time_set = pos(end,1) * time_set;
    pos_set = nan(0,3);
    
    for time_set_i = 1:size(time_set, 1)
        for PC_i = 1:length(spike_timestamps)
            spike_set_temp = spike_timestamps{PC_i}...
                (spike_timestamps{PC_i} >= time_set(time_set_i,1) & ...
                spike_timestamps{PC_i} <= time_set(time_set_i,2));
            spike_set_PC{PC_i,1} = [spike_set_PC{PC_i,1}; spike_set_temp];
        end
        
        pos_temp = pos(pos(:,1) >= time_set(time_set_i,1) & ...
            pos(:,1) <= time_set(time_set_i,2), :);
        pos_set = [pos_set; pos_temp]; %#ok<AGROW>
    end
end

function [training_set, testing_set] = setBin(set_i, bin_size)
    if set_i == 1
        training_set = [set_i*bin_size 1];
    elseif set_i*bin_size == 1
        training_set = [0 (set_i-1)*bin_size];
    else
        training_set = [[0 (set_i-1)*bin_size]; [set_i*bin_size 1]];
    end
    testing_set = [bin_size*(set_i-1) bin_size*set_i];
end

function [fit_loss_PC_mae, predict_loss_PC_mae] = decodeIntegration(spike_timestamps_PC, pos, p_map, ...
        bin_size, drawfig, session_ID, cell_type)
    % bin_size = 0.2;
    [fit_loss_PC_mae, predict_loss_PC_mae] = deal(nan(1,ceil(1/bin_size)));
    for set_i = 1:ceil(1/bin_size)
        [training_set, testing_set] = setBin(set_i, bin_size);
        
        % divide dataset;
        [spike_training_PC, pos_training_PC] = divideData(spike_timestamps_PC, pos, training_set);
        [spike_testing_PC, pos_testing_PC] = divideData(spike_timestamps_PC, pos, testing_set);
        spike_testing_PC(cellfun(@isempty, spike_training_PC(:))) = [];
        spike_training_PC(cellfun(@isempty, spike_training_PC(:))) = [];
        if isempty(spike_training_PC) || isempty(spike_testing_PC)
            [fit_loss_PC_mae(1,set_i), predict_loss_PC_mae(1,set_i)] = deal(nan);
            continue;
        end
        
        % decode;
        % PC;
        decoder_m = GNBDecoder('p_map', p_map, 't_size', 1, 't_smooth', 10);
        decoder_fit_PC = decoder_m.fit(pos_training_PC, spike_training_PC);
        decoder_predict_PC = decoder_m.predict(pos_testing_PC, spike_testing_PC);
        
        fit_loss_PC_mae(1,set_i) = decoder_fit_PC.loss.mae;
        predict_loss_PC_mae(1,set_i) = decoder_predict_PC.loss.mae;
        % fit_loss_PC{PC_type_i}(folder_i,:) = [decoder_fit_PC.loss.mae, decoder_fit_PC.loss.mse];
        % predict_loss_PC{PC_type_i}(folder_i,:) = [decoder_predict_PC.loss.mae, decoder_predict_PC.loss.mse];
        
        % figure;
        if drawfig
            plotDecodedTraj(decoder_fit_PC, decoder_predict_PC, session_ID, cell_type);
            saveas(gcf, ['G:\ZX\Data_temp\Results_MEC_temp\decode\', ...
                session_ID, '_', cell_type, '_', num2str(set_i), '.png']);
            close all;
        end
    end
end

function [fit_loss_PC_mae_shuffle, predict_loss_PC_mae_shuffle] = decodeShuffle(spike_timestamps_PC, pos, p_map, bin_size)
    shuffle_num = 100; % 1000;
    [fit_loss_PC_mae_shuffle, predict_loss_PC_mae_shuffle] = deal(cell(shuffle_num,1));
    
    % position shuffle;
    time_30s = find(pos(:,1) >= 30, 1);
    pos_shuffle_i = time_30s + randi((size(pos,1) - time_30s*2), shuffle_num,1);
    pos_shuffle = cell(shuffle_num,1);
    
    for n = 1:shuffle_num
        pos_shuffle{n} = [pos(pos_shuffle_i(n)+1:end, :);pos(1:pos_shuffle_i(n),:)];
        pos_shuffle{n}(:,1) = pos(:,1);
    end
    
    parfor n = 1:shuffle_num
        [fit_loss_PC_mae_shuffle{n,1}, predict_loss_PC_mae_shuffle{n,1}] = ...
            decodeIntegration(spike_timestamps_PC, pos_shuffle{n}, p_map, bin_size, 0, '', '');
    end
end

function plotDecodedTraj(decoder_fit_PC, decoder_predict_PC, session_ID, PC_type)
    figure('Position', [400 250 1200 600]);
    subplot('Position', [0.1 0.55 0.35 0.35]);
    hold on;
    plot(decoder_fit_PC.t_i, decoder_fit_PC.x_i(1,:), 'Color',[0.5 0.5 0.5], 'LineWidth',1.5);
    plot(decoder_fit_PC.t_i, decoder_fit_PC.x_hat(1,:), 'Color',[1 0 0 0.3], 'LineWidth',1.5);
    hold off;
    % xlabel('time (sec)');
    ylabel('x (cm)');
    title(['session: ', session_ID, ' - ', PC_type, ' - fit MAE: ', num2str(decoder_fit_PC.loss.mae)], 'Interpreter','none');
    
    subplot('Position', [0.1 0.1 0.35 0.35]);
    hold on;
    plot(decoder_fit_PC.t_i, decoder_fit_PC.x_i(2,:), 'Color',[0.5 0.5 0.5], 'LineWidth',1.5);
    plot(decoder_fit_PC.t_i, decoder_fit_PC.x_hat(2,:), 'Color',[1 0 0 0.3], 'LineWidth',1.5);
    hold off;
    xlabel('time (sec)');
    ylabel('y (cm)');
    
    
    subplot('Position', [0.55 0.55 0.35 0.35]);
    hold on;
    plot(decoder_predict_PC.t_i, decoder_predict_PC.x_i(1,:), 'Color',[0.5 0.5 0.5], 'LineWidth',1.5);
    plot(decoder_predict_PC.t_i, decoder_predict_PC.x_hat(1,:), 'Color',[1 0 0 0.3], 'LineWidth',1.5);
    hold off;
    % xlabel('time (sec)');
    ylabel('x (cm)');
    title([PC_type, ' - predict MAE: ', num2str(decoder_predict_PC.loss.mae)], 'Interpreter','none');
    
    subplot('Position', [0.55 0.1 0.35 0.35]);
    hold on;
    plot(decoder_predict_PC.t_i, decoder_predict_PC.x_i(2,:), 'Color',[0.5 0.5 0.5], 'LineWidth',1.5);
    plot(decoder_predict_PC.t_i, decoder_predict_PC.x_hat(2,:), 'Color',[1 0 0 0.3], 'LineWidth',1.5);
    hold off;
    xlabel('time (sec)');
    ylabel('y (cm)');
    legend(["Predicted", "Observed"]);
end

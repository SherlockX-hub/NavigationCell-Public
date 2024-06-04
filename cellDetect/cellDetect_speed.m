%% cellDetect_speed.m
% This script is used to find multiple speed cells and AHV cells.

% Created by Xiang Zhang, Nov., 2023.

clear;
tic;

%% code path;
addpath(genpath('G:\ZX\Codes\BNT-master'));

%% data;
dir_name = 'G:\ZX\Data_temp\49-20230108\49-20230108-4';
load([dir_name, '/NeuronActivity.mat'], 'NeuronActivity');
load([dir_name, '/behav.mat'], 'behav');

%% parameters;
do_cell_analysis = struct;
do_cell_analysis.speed_cell = 1;

cell_name = struct;
cell_name.speed_cell = 'speed_cell';

timestamps_min = min(NeuronActivity.timestamps(1), behav.timestamps(1));
calcium_time = seconds(NeuronActivity.timestamps - timestamps_min);
behav_time = seconds(behav.timestamps - timestamps_min);
calcium_event = NeuronActivity.Event_filtered_exp2;
pos_self = behav.position{1};
pos = [behav_time, pos_self];
hd_dir = behav.hdDir{1};

p = struct;
speed_sample_time = 0.02; p.speed2.speed_sample_time = speed_sample_time;
speed_win = 100; p.speed2.speed_win = speed_win; % time window;
p.speed2.speed_range = [0 0 0];
p.speed2.speed_bin = 1;

time_all = 0:speed_sample_time:min(calcium_time(end), behav_time(end));

% smooth the firing rate;
speed_firing_rate = NeuronActivity.Event_raw_exp2;
speed_firing_rate = interp1(calcium_time, speed_firing_rate, time_all'); % 20
p.speed2.speed_firing_rate = cell2mat(arrayfun(@(x) general.smoothGauss(speed_firing_rate(:,x), 0.4/speed_sample_time), ...
    1:size(speed_firing_rate,2), 'UniformOutput', false));

speed_filter = [1 2.5 0];
shuffle_num = 1000;
threshold_prc = 95;
fig_fmt = {'.fig', '.png'};

color_cell = struct;
color_cell.cell = [0.8 0 0];
color_cell.speed_cell = [];

% draw_fig = 1;
draw_fig.cell = 0;
SFP = NeuronActivity.SFP;

% main function;
%% original speed;
% calculate and smooth speed;
pos_smooth = interp1(behav_time, pos, time_all);
speed_x = [nan(speed_win/2,1); pos_smooth(speed_win+1:end,2) - pos_smooth(1:end-speed_win,2); nan(speed_win/2,1)];
speed_y = [nan(speed_win/2,1); pos_smooth(speed_win+1:end,3) - pos_smooth(1:end-speed_win,3); nan(speed_win/2,1)];
p.speed2.speed_input = sqrt(speed_x.^2 + speed_y.^2) / (speed_win * speed_sample_time);
p.speed2.speed_range = [0 0 0];
save_folder = [dir_name, '/speed_cell/speed_cell']; if ~exist(save_folder, 'dir'), mkdir(save_folder); end
cellDetect(do_cell_analysis, cell_name, calcium_time, calcium_event, pos, hd_dir, save_folder, p, ...
    'speed_filter',speed_filter, 'shuffle_num',shuffle_num, 'threshold_prc',threshold_prc, ...
    'draw_fig',draw_fig, 'SFP',SFP, 'fig_fmt',fig_fmt, 'color_cell',color_cell);

% reverse;
load([save_folder, '/speed_cell/speed_score.mat'], 'speed_score');
load([save_folder, '/speed_cell/speed_score_shuffle.mat'], 'speed_score_shuffle');
threshold_low = prctile(speed_score_shuffle(1:1000,:), 100-threshold_prc);
speed_cell_reverse = find((speed_score - threshold_low') < 0);
save([save_folder, '/speed_cell/speed_cell_reverse.mat'], 'speed_cell_reverse');

% angular head velocity;
%% absolute;
% ahv = analyses.angularHeadVelocity([behav_time, hd_dir]);
% p.speed2.speed_input = interp1(behav_time, ahv, time_all');

% speed_win = 100; p.speed2.speed_win = speed_win;
p.speed2.speed_bin = 3;
hd_dir_smooth = interp1(behav_time, hd_dir, time_all');
hd_dir_end = hd_dir_smooth(speed_win+1:end,1);
hd_dir_start = hd_dir_smooth(1:end-speed_win,1);
delta_hd_dir = [nan(speed_win/2,1); angleDiffer(hd_dir_end, hd_dir_start); nan(speed_win/2,1)];
p.speed2.speed_input = delta_hd_dir / (speed_win * speed_sample_time);
p.speed2.speed_range = [0 0 0];
save_folder = [dir_name, '/speed_cell/ahv_cell']; if ~exist(save_folder, 'dir'), mkdir(save_folder); end
cellDetect(do_cell_analysis, cell_name, calcium_time, calcium_event, pos, hd_dir, save_folder, p, ...
    'speed_filter',speed_filter, 'shuffle_num',shuffle_num, 'threshold_prc',threshold_prc, ...
    'draw_fig',draw_fig, 'SFP',SFP, 'fig_fmt',fig_fmt, 'color_cell',color_cell);

% reverse;
load([save_folder, '/speed_cell/speed_score.mat'], 'speed_score');
load([save_folder, '/speed_cell/speed_score_shuffle.mat'], 'speed_score_shuffle');
threshold_low = prctile(speed_score_shuffle(1:1000,:), 100-threshold_prc);
speed_cell_reverse = find((speed_score - threshold_low') < 0);
save([save_folder, '/speed_cell/speed_cell_reverse.mat'], 'speed_cell_reverse');

%% cw or ccw avh;
% ccw;
hd_dir_sign = sign(cell2mat(arrayfun(@(a,b,c,d) det([a,b; c,d]), ...
    cosd(hd_dir_end), sind(hd_dir_end), cosd(hd_dir_start), sind(hd_dir_start), 'UniformOutput', false)));
delta_hd_dir = [nan(speed_win/2,1); hd_dir_sign .* angleDiffer(hd_dir_end, hd_dir_start); nan(speed_win/2,1)];
speed_input = delta_hd_dir / (speed_win * speed_sample_time);
p.speed2.speed_input = speed_input;
p.speed2.speed_input(speed_input < 0) = nan;
p.speed2.speed_range = [0 0 0];
save_folder = [dir_name, '/speed_cell/ahv_cell_ccw']; if ~exist(save_folder, 'dir'), mkdir(save_folder); end
cellDetect(do_cell_analysis, cell_name, calcium_time, calcium_event, pos, hd_dir, save_folder, p, ...
    'speed_filter',speed_filter, 'shuffle_num',shuffle_num, 'threshold_prc',threshold_prc, ...
    'draw_fig',draw_fig, 'SFP',SFP, 'fig_fmt',fig_fmt, 'color_cell',color_cell);

% reverse;
load([save_folder, '/speed_cell/speed_score.mat'], 'speed_score');
load([save_folder, '/speed_cell/speed_score_shuffle.mat'], 'speed_score_shuffle');
threshold_low = prctile(speed_score_shuffle(1:1000,:), 100-threshold_prc);
speed_cell_reverse = find((speed_score - threshold_low') < 0);
save([save_folder, '/speed_cell/speed_cell_reverse.mat'], 'speed_cell_reverse');

%% cw;
p.speed2.speed_input = speed_input;
p.speed2.speed_input(speed_input > 0) = nan;
p.speed2.speed_input = abs(p.speed2.speed_input);
p.speed2.speed_range = [0 0 0];
save_folder = [dir_name, '/speed_cell/ahv_cell_cw']; if ~exist(save_folder, 'dir'), mkdir(save_folder); end
cellDetect(do_cell_analysis, cell_name, calcium_time, calcium_event, pos, hd_dir, save_folder, p, ...
    'speed_filter',speed_filter, 'shuffle_num',shuffle_num, 'threshold_prc',threshold_prc, ...
    'draw_fig',draw_fig, 'SFP',SFP, 'fig_fmt',fig_fmt, 'color_cell',color_cell);

% reverse;
load([save_folder, '/speed_cell/speed_score.mat'], 'speed_score');
load([save_folder, '/speed_cell/speed_score_shuffle.mat'], 'speed_score_shuffle');
threshold_low = prctile(speed_score_shuffle(1:1000,:), 100-threshold_prc);
speed_cell_reverse = find((speed_score - threshold_low') < 0);
save([save_folder, '/speed_cell/speed_cell_reverse.mat'], 'speed_cell_reverse');

%% all ahv range;
p.speed2.speed_input = speed_input;
p.speed2.speed_range = [1 -180 180];
save_folder = [dir_name, '/speed_cell/ahv_cell_ccw_cw']; if ~exist(save_folder, 'dir'), mkdir(save_folder); end
cellDetect(do_cell_analysis, cell_name, calcium_time, calcium_event, pos, hd_dir, save_folder, p, ...
    'speed_filter',speed_filter, 'shuffle_num',shuffle_num, 'threshold_prc',threshold_prc, ...
    'draw_fig',draw_fig, 'SFP',SFP, 'fig_fmt',fig_fmt, 'color_cell',color_cell);

% reverse;
load([save_folder, '/speed_cell/speed_score.mat'], 'speed_score');
load([save_folder, '/speed_cell/speed_score_shuffle.mat'], 'speed_score_shuffle');
threshold_low = prctile(speed_score_shuffle(1:1000,:), 100-threshold_prc);
speed_cell_reverse = find((speed_score - threshold_low') < 0);
save([save_folder, '/speed_cell/speed_cell_reverse.mat'], 'speed_cell_reverse');

toc;
%% cellDetect_run.m
% This script is an example of how function cellDetect.m works.
% Some parameters are optional and variable based on users' need.
%
% Created by Xiang Zhang, March, 2023.

clear;
tic;

%% code path;
addpath(genpath('G:\ZX\Codes\BNT-master'));

%% data;
dir_name = pwd;
load([dir_name, '/My_V4_Miniscope/ms.mat'], 'ms');
load([dir_name, '/My_V4_Miniscope/SFP.mat'], 'SFP');
load([dir_name, '/My_WebCam/behav.mat'], 'behav');

%% timestamps correct;

%% parameters;
do_cell_analysis = struct;
do_cell_analysis.place_cell = 1;
do_cell_analysis.grid_cell = 1;
do_cell_analysis.border_cell = 1;
do_cell_analysis.head_direction_cell = 1;
do_cell_analysis.speed_cell = 0;

cell_name = struct;
cell_name.place_cell = 'place_cell';
cell_name.grid_cell = 'grid_cell';
cell_name.border_cell = 'border_cell';
cell_name.head_direction_cell = 'head_direction_cell';
cell_name.speed_cell = 'speed_cell';

calcium_time = ms.time/1000;
behav_time = behav.time/1000;
calcium_event = ms.DeconvSignals > repmat(3*std(ms.DeconvSignals), size(ms.DeconvSignals,1),1);
pos = [behav_time, behav.position_smooth{1}];
% hd_dir = behav.hdDir_smooth{1};
hd_dir_eul = rad2deg(quat2eul(ms.HeadOrientation));
hd_dir_ms = mod(-hd_dir_eul(:,1), 360);
hd_dir = mod(rad2deg(angle(interp1(ms.time, exp(deg2rad(hd_dir_ms)*1i), behav.time))),360);
save_folder = [dir_name, '/self_cell'];

% minimum event number;
min_event_num = 100;

% spatial dimension;
spatial_dimension = 2;

p = struct;
% map parameters;
p.map.datatime = 's'; % 'msec', 's'
p.map.binWidth = 2;
p.map.smooth = 2;
p.map.minTime = 0;
p.map.maxGap = 0.3;
p.map.limits = [0 behav.trackLength(1) 0 behav.trackLength(2)];
p.map.blanks = 'on';

% field parameters;
p.field.binWidth = 2;
p.field.minPeak = 0.5;
p.field.threshold = 0.6;

% grid cell parameters;
% p.gird. = ;
p.grid2.radii = [5 5];

% head direction cell parameters;
p.hd.binWidth = 3;
p.hd.smooth = 2;
p.hd2.sampleTime = mean(diff(behav_time));
p.hd2.percentile = 50;
% p.hd2.trajectoryNorm = 1/4;
p.hd2.trajectoryPlot = {'Color', [0.28 0.6 0.75294]};

% border cell parameters;
% p.border. = ;

%{
% speed cell preprocession;
% p.speed2.span = 20;
speed_sample_time = 0.02; p.speed2.speed_sample_time = speed_sample_time;
speed_win = 100; p.speed2.speed_win = speed_win; % time window;
p.speed2.speed_range = [0 0 0];
p.speed2.speed_bin = 1;
% generate the same time line for calcium event and speed;
time_all = 0:speed_sample_time:min(calcium_time(end), behav_time(end));

% smooth the firing rate;
speed_firing_rate = ms.DeconvSignals;
speed_firing_rate = interp1(calcium_time, speed_firing_rate, time_all); % 20
p.speed2.speed_firing_rate = cell2mat(arrayfun(@(x) general.smoothGauss(speed_firing_rate(:,x), 0.4/speed_sample_time), ...
    1:size(speed_firing_rate,2), 'UniformOutput', false));

% calculate and smooth speed;
pos_smooth = interp1(behav_time, pos, time_all);
speed_x = [nan(speed_win/2,1); pos_smooth(speed_win+1:end,2) - pos_smooth(1:end-speed_win,2); nan(speed_win/2,1)];
speed_y = [nan(speed_win/2,1); pos_smooth(speed_win+1:end,3) - pos_smooth(1:end-speed_win,3); nan(speed_win/2,1)];
p.speed2.speed_input = sqrt(speed_x.^2 + speed_y.^2) / (speed_win * speed_sample_time);
%}

% map speed filter;
speed_filter = [1 8 0];

% shuffle;
shuffle_num = 1000;
threshold_prc = 95; % threshold_prc.head_direction_cell = 99;

% figure;
draw_fig = 1; % draw_fig.cell = 0;
p.behav_limit = [0 behav.trackLength(1) 0 behav.trackLength(2)]; % p.map.limits;

% load([dir_name, '/Fall.mat'], 'stat', 'ops', 'iscell');
% SFP = suite2pSFP(stat, ops, iscell);
fig_fmt = {'.fig', '.png'};

color_cell = struct;
color_cell.cell = [0.8 0 0];
color_cell.place_cell = [];
color_cell.grid_cell = [];
color_cell.border_cell = [];
color_cell.head_direction_cell = [];
color_cell.speed_cell = [];

%% main function;
cellDetect(do_cell_analysis, cell_name, calcium_time, calcium_event, pos, hd_dir, save_folder, p, ...
    'min_event_num',min_event_num, 'spatial_dimension',spatial_dimension, ...
    'speed_filter',speed_filter, 'shuffle_num',shuffle_num, 'threshold_prc',threshold_prc, ...
    'draw_fig',draw_fig, 'SFP',SFP, 'fig_fmt',fig_fmt, 'color_cell',color_cell);

toc;
%% GNBDecoder.m
% This is a Gaussian Naive Bayes decoder based on Zhang KC et al., 1998
% Modified from Kai Gao.

% Input:
%        pos: a p*3 matrix, [time x y];
%             or a p*2 matrix, [time x] in one-dimension;
%        spike_time_stamps: a cell, containing spike timestamps in seconds of all cells;
%
% Optional input:
%        p_map: a struct, parameters of rate map;
%               (details seen in map.m of BNT package)
%        t_size: a scale, width of time bins, (default, 1);
%        t_smooth: a scale, smooth size in time bins (default = 10);
%
% Output:
%        estimate: a struct, containing decoding results;
%                  estimate.t_i: a vector, mean time of each time bin;
%                  estimate.x_i: a vector, mean position at each time bin;
%                  estimate.x_hat: a vector, predicted position at each time bin;
%                  estimate.loss: a struct, m.a.e. and m.s.e. of the result;
%
% Usage:
%        initiation: decoder_GNB = GNBDecoder(varargin);
%        training: result_training = decoder_GNB.fit(pos_training, spike_time_stamps_training);
%        testing: result_testing = decoder_GNB.predict(pos_testing, spike_time_stamps_testing);

% Modified by Xiang Zhang, Sept., 2023.

% Modification:
% change the input of parameters;
% merge 1D and 2D version;
% correct the calculation in mae;

classdef GNBDecoder < handle
    
    properties
        p_map;
        t_size;
        t_smooth;
        numNeurons;
        f_i; % firing rate;
        p_x_mask;
        id2position; % the postion value in each bin of rate map;
    end
    
    methods
        function obj = GNBDecoder(varargin)
            inp = inputParser;
            addParameter(inp, 'p_map', struct);
            addParameter(inp, 't_size', 1);
            addParameter(inp, 't_smooth', 1);
            
            parse(inp, varargin{:});
            p_map = inp.Results.p_map; obj.p_map = collectPara(p_map);
            obj.t_size = inp.Results.t_size;
            obj.t_smooth = inp.Results.t_smooth;
            
            % check dependent packages
            if ~exist('BNT', 'dir') && ~exist('BNT-master', 'dir'), error('Please add BNT package to Path first.'); end
        end
        
        function estimate = fit(obj, pos, spike_time_stamps)
            if size(pos,2) < 2 || size(pos,2) > 3, error('Wrong input of ''pos''.'); end
            if ~iscell(spike_time_stamps), error('Wrong input of ''spike_time_stamps''.'); end
            if any(cellfun(@isempty, spike_time_stamps)), error(['There is a cell whose firing rate is zero in training set. ', ...
                    'Decoding terminated. Please double check and remove it.']); end
            
            obj.numNeurons = length(spike_time_stamps);
            
            % binned time and pos;
            [estimate.t_i, estimate.x_i, bin_idx, bin_number] = obj.discretize_traj(pos);
            
            % estimate P(x) for the trial, and f_i(x) for each cell
            obj.f_i = cell(obj.numNeurons, 1);
            n_i = zeros(obj.numNeurons, bin_number);
            h = waitbar(0, 'training on neurons...');
            for cell_i = 1:obj.numNeurons
                waitbar(cell_i/obj.numNeurons, h, ['training on neurons...' num2str(cell_i) '/' num2str(obj.numNeurons)]);
                % ------------- training set --------------
                r_timestamp = pos(knnsearch(pos(:, 1), spike_time_stamps{cell_i}(:)), 1);

                % f_i(x)
                map = analyses.map(pos, r_timestamp, obj.p_map{:});
                map.z(isnan(map.z)) = 0;
                obj.f_i{cell_i} = map.z(:);

                % n_i(t)
                tn_i = arrayfun(@(x) sum(pos(bin_idx == x, 1)' == r_timestamp, 'all'), 1:bin_number);
                n_i(cell_i, :) = smoothdata(tn_i, 'gaussian', obj.t_smooth);
            end
            close(h);
            
            % p(x)
            map.time(isnan(map.time)) = 0;
            obj.p_x_mask = map.time(:); obj.p_x_mask(map.time~=0) = 1; obj.p_x_mask = log(obj.p_x_mask);
            % p_x_prior = rate_map.time(:) / sum(rate_map.time(:));
            x = map.x(1:(end-1)) + (map.x(2)-map.x(1)) / 2;
            if ~isempty(map.y)
                y = map.y(1:(end-1)) + (map.y(2)-map.y(1)) / 2;
                [pos_x, pos_y] = meshgrid(x, y);
                obj.id2position = arrayfun(@(x,y) [x,y], pos_x, pos_y, 'UniformOutput', false);
            else
                obj.id2position = arrayfun(@(x) x, x, 'UniformOutput', false);
            end
            
            [estimate.x_hat, estimate.loss] = decode(obj, estimate.x_i, obj.p_x_mask, n_i, 'decoding on training set...');
        end
        
        function estimate = predict(obj, pos, spike_time_stamps)
            if size(pos,2) < 2 || size(pos,2) > 3, error('Wrong input of ''pos''.'); end
            if ~iscell(spike_time_stamps), error('Wrong input of ''spike_time_stamps''.'); end
            if obj.numNeurons ~= length(spike_time_stamps)
                error('Neuron number in test set doesn''t match that in training set.'); end
            
            % split test set time bin
            [estimate.t_i, estimate.x_i, bin_idx, bin_number] = obj.discretize_traj(pos);
            
            % estimate P(x) for the trial, and f_i(x) for each cell
            n_i = zeros(obj.numNeurons, bin_number);
            h = waitbar(0, 'preprocessing test set...');
            for cell_i = 1:obj.numNeurons
                waitbar(cell_i/obj.numNeurons, h, ['preprocessing test set...' num2str(cell_i) '/' num2str(obj.numNeurons)]);
                % ------------- test set --------------
                r_timestamp = pos(knnsearch(pos(:, 1), spike_time_stamps{cell_i}(:)), 1);
                
                % n_i(t)
                cn_i = arrayfun(@(x) sum(pos(bin_idx == x, 1)' == r_timestamp, 'all'), 1:bin_number);
                n_i(cell_i, :) = smoothdata(cn_i, 'gaussian', obj.t_smooth);
            end
            close(h);
            
            [estimate.x_hat, estimate.loss] = decode(obj, estimate.x_i, obj.p_x_mask, n_i, 'decoding on test set...');
        end
        
        function [t_i, x_i, bin_idx, bin_number] = discretize_traj(obj, pos)
            t_min = nanmin(pos(:,1)); t_min = t_min - mod(t_min, obj.t_size);
            t_max = nanmax(pos(:,1)); t_max = t_max - mod(t_max, obj.t_size) + obj.t_size;
            [~, ~, bin_idx] = histcounts(pos(:,1), t_min:obj.t_size:t_max);
            
            bin_number = max(bin_idx);
            t_i = arrayfun(@(x) nanmean(pos(bin_idx == x, 1), 1), (1:bin_number)'); % mean time within a time bin;
            x_i = cell2mat(arrayfun(@(x) nanmean(pos(bin_idx == x, 2:end), 1), (1:bin_number)', ...
                'UniformOutput',false)); % mean position within a time bin;
        end
        
        function [x_hat, loss] = decode(obj, x_i, p_x_mask, n_i, hint)
            bin_number = size(n_i, 2);
            x_hat = nan(bin_number, length(obj.id2position{1}));
            h = waitbar(0, hint);
            for bin_id = 1:bin_number
                waitbar(bin_id/bin_number, h, [hint num2str(bin_id) '/' num2str(bin_number)])

                p_x_posterior = zeros(length(p_x_mask), 1);
                cn_i = n_i(:, bin_id);
                tf_i = obj.f_i;
                for cell_id = 1:obj.numNeurons
                    % assume Poisson neurons
                    p_x = p_x_mask + cn_i(cell_id).*log(tf_i{cell_id}+1e-50) - tf_i{cell_id}*obj.t_size;  % assume prior == 1
                    % p_x = log(p_x_prior .* (tf_i{cell_id}).^cn_i(cell_id).*exp(-tf_i{cell_id}*obj.t_size));
                    p_x_posterior = p_x_posterior + p_x;
                end

                [~, decoded_pos_id] = max(p_x_posterior);
                x_hat(bin_id, :) = obj.id2position{decoded_pos_id};
            end
            close(h)
            
            loss.mae = nanmean(vecnorm(x_i - x_hat, 2, 2));
            loss.mse = nanmean(vecnorm(x_i - x_hat, 2, 2).^2);
        end
        
    end
end

%% functions;
function para = collectPara(p)
    para = {};
    if isempty(p), return; end
    p_var = fieldnames(p); if isempty(p_var), return; end
    for var_i = 1:length(p_var)
        para{1,end+1} = p_var{var_i}; %#ok<AGROW>
        para{1,end+1} = p.(p_var{var_i}); %#ok<AGROW>
    end
end

function filename = genealogy_experiment(n_reps, job_id, solver, dir_id)
% GENEALOGY_EXPERIMENT Run the genealogy experiment to produce data for Fig. 5.
%
% Input:
% 	n_reps -- Number of times the experiment is run with a random initialization.
% 	job_id -- Job ID number appended to the name of the output directory
% 		  (set it to, e.g., 0 unless you run a batch job).
% 	solver -- Should be either 0 (for Natalie) or 1 (for NetAlignMP++).
% 	dir_id -- ID number for the experiment appended to the output directory.
%
% Author: Eric Malmi (eric.malmi@gmail.com)

addpath('netalign/matlab');

if nargin < 4
    dir_id = 0;
end
if nargin < 3
    solver = 1;
end
if nargin < 2
    job_id = 0;
end
if nargin < 1
    n_reps = 1;
end

% Set random seed.
rng(now() + job_id);

subdir = strcat('experiment_results/genealogy_run', int2str(dir_id));
mkdir(subdir);
filename = strcat(subdir, '/res', int2str(job_id), '_solver', ...
    int2str(solver), '.mat');

% Problem parameters.
data_dir = 'genealogy_data/';%_small30/';
config = load(strcat(data_dir, 'config.mat'));
n = double(config.num_source_nodes);
m = double(config.num_target_nodes);
data_getter = @(repetition_idx) get_genealogy_problem( ...
    strcat(data_dir, 'problem', int2str(job_id+repetition_idx-1)), n, m);

% Solver parameters.
k = 30;
a = 1;
b = 1;
gamma = 0.1;
stepm = 20;
dtype = 2;
verbose = false;
max_iters = 300;

% Experiment parameters
method_names = {'TopMatchings', 'LCCL', 'Margin', ... %'Degree', 'MinDegree', ...
                'Betweenness', 'Random', 'TopMatchings10Batch', ...
                'Margin10Batch'};
query_counts = cell(length(method_names), 1);
for i = 1:length(method_names)
    if strfind(method_names{i}, 'Batch')
        query_counts{i} = 0:10:990;
    else
        query_counts{i} = 0:990;
    end
end

[accs, ts, query_ts, cum_accs, align_ts] = run_experiment( ...
    method_names, query_counts, n_reps, n, data_getter, solver, a, b, ...
    gamma, stepm, dtype, max_iters, verbose, k);
mean_query_ts = mean(query_ts, 3)
mean_ts = mean(ts, 3)
mean_cum_accs = mean(cum_accs, 3)
mean_accs = mean(accs, 3)

save(filename);
fprintf('Wrote workspace to: %s\n', filename);

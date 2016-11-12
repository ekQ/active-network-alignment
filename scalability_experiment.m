function filename = scalability_experiment(n_reps, job_id, solver, dir_id)
% SCALABILITY_EXPERIMENT Run the scalability experiment (Table 1).
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
rng(mod(now()*1e5,1e5) + job_id);

subdir = strcat('experiment_results/scalability_run', int2str(dir_id));
mkdir(subdir);
filename = strcat(subdir, '/res', int2str(job_id), '_solver', ...
    int2str(solver), '.mat');

% Problem parameters.
p_keep_edge = 0.4;
density_multiplier = 1.5;
n_duplicates = 30;

% Solver parameters.
k = 30;
a = 1;
b = 1;
gamma = 0.1;
stepm = 20;
dtype = 2;
verbose = true;
max_iters = 300;

% Experiment parameters
method_names = {'TopMatchings'};
query_counts = cell(length(method_names), 1);
for i = 1:length(method_names)
    query_counts{i} = [0, 1];
end

ns = [100, 1000, 10000];
all_query_ts = zeros(length(ns), n_reps);
all_align_ts = zeros(length(ns), n_reps);
for i = 1:length(ns)
    n = ns(i);
    data_getter = @(repetion_idx) get_synthetic_problem(n, 2, p_keep_edge, ...
        density_multiplier, n_duplicates);
    [accs, ts, query_ts, cum_accs, align_ts] = run_experiment( ...
        method_names, query_counts, n_reps, n, data_getter, solver, a, b, ...
        gamma, stepm, dtype, max_iters, verbose, k);
    mean_query_ts = mean(query_ts, 3)
    mean_align_ts = mean(align_ts, 3)
    mean_accs = mean(accs, 3)
    all_query_ts(i,:) = query_ts(1,2,:);
    all_align_ts(i,:) = align_ts(1,1,:);
end

qt = median(all_query_ts, 2)
at = median(all_align_ts, 2)

save(filename);
fprintf('Wrote workspace to: %s\n', filename);

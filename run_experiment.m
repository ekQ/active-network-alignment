function [accs, ts, query_ts, cum_accs, align_ts] = run_experiment( ...
        method_names, all_query_counts, n_reps, n, data_getter, solver, ...
        a, b, gamma, stepm, dtype, max_iters, verbose, k)
% RUN_EXPERIMENT A generic experiment template for the experiments of the
% paper.
%
% Author: Eric Malmi (eric.malmi@gmail.com)

max_query_count_len = -1;
assert(length(all_query_counts) == length(method_names));
for i = 1:length(all_query_counts)
    len = length(all_query_counts{i});
    if len > max_query_count_len
        max_query_count_len = len;
    end
end
accs = zeros(length(method_names), max_query_count_len, n_reps);
cum_accs = zeros(length(method_names), max_query_count_len, n_reps);
ts = zeros(length(method_names), max_query_count_len, n_reps);
query_ts = zeros(length(method_names), max_query_count_len, n_reps);
align_ts = zeros(length(method_names), max_query_count_len, n_reps);
for r = 1:n_reps
    [A, B, L0, groundtruth] = data_getter(r);
    accs_slice = zeros(length(method_names), max_query_count_len);
    for i = 1:length(method_names)
        accs_slice_slice = zeros(1, max_query_count_len);
        method = method_names{i};
        query_counts = all_query_counts{i};
        n_queried = 0;
        already_queried = logical(zeros(n,1));
        for j = 1:length(query_counts)
            t0 = tic;
            L = L0;
            qc = query_counts(j);
            % Query if needed
            if n_queried < qc
                n_query_next = qc - n_queried;
                next_nodes = query_nodes(A, bestL, n_query_next, method, ...
                                         already_queried, k, verbose);
                already_queried(next_nodes) = 1;
                n_queried = n_queried + n_query_next;
                % Set fixed node rows to zero except for the true match and
                % the same for the columns of the true matches.
                true_idxs = sub2ind(size(L), find(already_queried), ...
                                    groundtruth(already_queried));
                L(already_queried, :) = 0;
                L(:, groundtruth(already_queried)) = 0;
                L(true_idxs) = 1;
            end
            query_ts(i, j, r) = toc(t0);
            % Ugly hack.
            if length(query_counts) == 2 && j == 2 && qc == 1
                continue
            end
            t0_align = tic;
            [S,w,li,lj] = netalign_setup(A,B,L);
            if solver == 0
                [~, ~, ~, ma, mb, prob] = netalignmr(S, w, a, b, li, lj, ...
                        gamma, stepm, 1, max_iters, verbose);
            elseif solver == 1
                [~, ~, ~, ~, ma, mb, prob] = netalignscbp(S, w, a, b, li, ...
                        lj, gamma, dtype, max_iters, verbose);
            end
            align_ts(i, j, r) = toc(t0_align);
            bestL = csr_to_sparse(prob.rp, prob.ci, prob.ai);
            bestL(:,size(B,1)+1:end) = [];
            % Compute accuracy on unqueried nodes.
            pred = -1 * ones(n,1);
            pred(ma) = mb;
            idxs = find(~already_queried);
            acc = sum(groundtruth(idxs) == pred(idxs)) / length(idxs);
            accs_slice_slice(j) = acc;

            % Compute accuracy on all nodes.
            pred(already_queried) = groundtruth(already_queried);  % Just in case.
            cum_accs(i, j, r) = sum(groundtruth == pred) / length(pred);
            ts(i, j, r) = toc(t0);
            fprintf('rep=%d, method=%s, qc=%d, acc=%.4f.\n', ...
                    r, method, qc, acc);
        end
        accs_slice(i, :) = accs_slice_slice;
    end
    accs(:, :, r) = accs_slice;
end

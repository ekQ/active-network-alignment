function [A, B, L, matches] = get_synthetic_problem(n, m, p_keep_edge, ...
    density_multiplier, n_duplicates)
% GET_SYNTHETIC_PROBLEM Generate a synthetic problem instance.
%
% A: adjacency matrix of the source graph.
% B: adjacency matrix of the target graph.
% L: candidate matches and their unary similarities.
% matches: groundtruth

edges = preferential_attachment_opt(n, m);
A = sparse(n,n);
A(sub2ind(size(A), edges(:,1), edges(:,2))) = 1;
B = A;

density = nnz(B) / numel(B);
p_add_edge = density * (density_multiplier - 1);

fast_method = true;

% Let's remove some edges.
if ~fast_method
    kept = triu(rand(n));
    kept = kept + kept';
    kept(kept < p_keep_edge) = 1;
    kept(kept < 1) = 0;
    B = B & kept;
else
    B = triu(B, 1);
    nnz_idxs = find(B);
    nnz_count = length(nnz_idxs);
    remove_count = floor((1 - p_keep_edge) * nnz_count);
    if remove_count >= 1
        remove_idxs = randsample(nnz_idxs, remove_count, false);
        B(remove_idxs) = 0;
    end
    B = B + B';
end
% Add edges.
if ~fast_method
    add = triu(rand(n));
    add = add + add';
    add(add < p_add_edge) = 1;
    add(add < 1) = 0;
    B = B | add;
    B = B - diag(diag(B));
else
    % Multiplier two is there since we only take the upper triangle in the end.
    add_count = floor(2 * nnz_count * (density_multiplier - 1));
    if add_count >= 1
        % For large graphs (n>10k) sample with replacement for efficiency).
        add_idxs = randsample(n*n, add_count, n > 10000);
        B(add_idxs) = 1;
        B = triu(B, 1);
        B = B + B';
    end
end

% Shuffle the order of nodes.
order = randperm(n);
[~, matches] = sort(order);
matches = matches';
B = B(order,:);
B = B(:,order);

n_labels = round(n / n_duplicates);
labelsA = randi(n_labels, n, 1);
labelsB = labelsA(order);

% Create the similarity matrix.
L = sparse(n, n);
for i = 1:n_labels
    idxsA = labelsA == i;
    idxsB = labelsB == i;
    L(idxsA, idxsB) = 1;
end

%vals = sort(unique(labelsA));
%count = histc(labelsA, vals);
%hist(count,unique(count))

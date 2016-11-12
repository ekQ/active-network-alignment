function nodes = query_nodes(A, L, n_nodes, method, already_queried, k, ...
                             verbose)
% QUERY_NODES Return k nodes to be labeled by an oracle.
%
% nodes = query_nodes(A, L, k, method, already_queries) return k nodes to
% query based on the adjacency matrix of the source graph A, matching
% scores L, query method, and an indicator vector of already queried nodes.
%
% Author: Eric Malmi (eric.malmi@gmail.com)
if nargin < 6
    verbose = true;
end

n = size(A,1);
if strcmp(method, 'Degree')
    degrees = sum(A,2);
    certainties = -degrees;
elseif strcmp(method, 'MinDegree')
    degrees = sum(A,2);
    certainties = degrees;
elseif strcmp(method, 'Random')
    certainties = rand(n,1);
elseif strfind(method, 'Margin')
    [max_val, max_idx] = max(L, [], 2);
    L(sub2ind(size(L), 1:n, max_idx')) = -1e10;
    max_val2 = max(L, [], 2);
    certainties = max_val - max_val2;
elseif strcmp(method, 'Betweenness')
    G = graph(A);
    centralities = centrality(G, 'betweenness');
    certainties = -centralities;
elseif strfind(method, 'LCCL')
    [~, edges] = k_best_matchings(L, 1, verbose);
    matched_idxs = edges{1}(:,1);
    idxs = sub2ind(size(L), matched_idxs, edges{1}(:,2));
    match_scores = inf * ones(n,1);
    match_scores(matched_idxs) = L(idxs);
    certainties = match_scores;
elseif strfind(method, 'TopMatchings')
    [~, edges] = k_best_matchings(L, k, verbose);
    match_counts = zeros(n,size(L,2)+1);  % Last column denotes unmatched.
    match_counts(:,end) = k;
    for i = 1:k
        matched_idxs = edges{i}(:,1);
        idxs = sub2ind(size(match_counts), matched_idxs, edges{i}(:,2));
        match_counts(idxs) = match_counts(idxs) + 1;
        idxs2 = sub2ind(size(match_counts), matched_idxs, ...
            (size(L,2)+1) * ones(length(matched_idxs),1));
        match_counts(idxs2) = match_counts(idxs2) - 1;
    end
    % All row sums should be equal.
    row_sums = sum(match_counts,2);
    assert(min(row_sums) == max(row_sums));
    certainties = max(match_counts,[],2);
end
% Already queried are certain.
certainties(already_queried) = inf;
[~, order] = sort(certainties, 'ascend');
nodes = order(1:n_nodes);

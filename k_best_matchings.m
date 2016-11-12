function [solution_vals, solution_edges] = k_best_matchings(X, k, verbose)
% K_BEST_MATCHING Murty's top-k assignments algorithm.
%
% Author: Eric Malmi (eric.malmi@gmail.com)

if nargin < 3
    verbose = true;
end

assert(max(size(X,1), size(X,2)) < 2^16)

% Remove all-zero columns to speed-up matching computation.
active_columns = find(sum(X,1))';
X = X(:, active_columns);

nass = 0;
solution_vals = zeros(1,k);
solution_edges = cell(1,k);

nrows = size(X,1);
ncols = size(X,2);
[val, m1, m2] = bipartite_matching(X);
nodes = struct('active_rows', logical(ones(nrows,1)), ...
               'active_cols', logical(ones(ncols,1)), ...
               'fixed_edges', zeros(0,2,'uint16'), 'solved_edges', [m1, m2], ...
               'val', val, 'forbidden_edges', []);
for i = 1:k
    % Find node with the max val.
    max_val = -inf;
    max_idx = -1;
    for j = 1:length(nodes)
        if nodes(j).val > max_val
            max_val = nodes(j).val;
            max_idx = j;
        end
    end
    % Record a new solution.
    solution_vals(i) = max_val;
    max_node = nodes(max_idx);
    solution_edges{i} = [max_node.fixed_edges; max_node.solved_edges];
    solution_edges{i} = [solution_edges{i}(:,1) active_columns(solution_edges{i}(:,2))];
    % Split the max node unless we've already found k matchings.
    if i < k
        nodes_len = length(nodes);
        % Initialize the struct array.
        nodes(nodes_len+size(max_node.solved_edges, 1)).val = -1;
        for j = 1:size(max_node.solved_edges, 1)
            new_edge = max_node.solved_edges(j,:);
            new_forbidden_edges = max_node.forbidden_edges;
            new_forbidden_edges(end+1) = sub2ind(size(X), new_edge(1), ...
                                                 new_edge(2));
            edges_to_add = max_node.solved_edges(1:j-1,:);
            new_fixed_edges = [max_node.fixed_edges; edges_to_add];
            new_active_rows = max_node.active_rows;
            new_active_cols = max_node.active_cols;
            % Remove rows and columns corresponding to the new edges.
            new_active_rows(edges_to_add(:,1)) = 0;
            new_active_cols(edges_to_add(:,2)) = 0;
            % Construct the weight matrix.
            X2 = X;
            X2(new_forbidden_edges) = 0;
            X2 = X2(new_active_rows, new_active_cols);

            [new_val, m1, m2] = bipartite_matching(X2);

            % Add the cost of the fixed edges to new_val.
            new_val = new_val + ...
                sum(X(sub2ind(size(X), new_fixed_edges(:,1), ...
                              new_fixed_edges(:,2))));

            nass = nass + 1;

            % Create the new node.
            nodes(nodes_len+j).active_rows = new_active_rows;
            nodes(nodes_len+j).active_cols = new_active_cols;
            nodes(nodes_len+j).fixed_edges = new_fixed_edges;
            nodes(nodes_len+j).val = new_val;
            nodes(nodes_len+j).forbidden_edges = new_forbidden_edges;
            % Get matching edges.
            row_idxs2 = find(new_active_rows);
            col_idxs2 = find(new_active_cols);
            nodes(nodes_len+j).solved_edges = uint16([row_idxs2(m1), ...
                                                      col_idxs2(m2)]);
        end
    end
    if verbose
        fprintf('k=%d, num_nodes=%d\n', i, length(nodes));
    end
    % Remove the max node.
    nodes(max_idx) = [];
end

function [A, B, L, matches] = get_social_network_problem(n_duplicates)
% GET_SOCIAL_NETWORK_PROBLEM Read a social network alignment problem and
% generate the labels of the nodes.
%
% A: adjacency matrix of the source graph.
% B: adjacency matrix of the target graph.
% L: candidate matches and their unary similarities.
% matches: groundtruth

E = load('social_network_data/CS-Aarhus_multiplex.edges');
% Just take the first two layers.
E(E(:,1) > 2, :) = [];
% Number of nodes is equal to the maximum node ID.
n = max(E(:));

% Construct adjacency matrices.
E_A = E(E(:,1)==2, 2:3);  % Facebook network.
E_B = E(E(:,1)==1, 2:3);  % Lunch network.
A = zeros(n);
B = zeros(n);
A(sub2ind([n, n], [E_A(:,1); E_A(:,2)], [E_A(:,2); E_A(:,1)])) = 1;
B(sub2ind([n, n], [E_B(:,1); E_B(:,2)], [E_B(:,2); E_B(:,1)])) = 1;

% Candidate matches matrix.
n_uniq_names = ceil(n / n_duplicates);
labels = randi(n_uniq_names, n, 1);
L = sparse(repmat(labels, 1, n) == repmat(labels', n, 1));

% Shuffle the order of nodes.
order = randperm(n);
[~, matches] = sort(order);
matches = matches';
B = B(order,:);
B = B(:,order);
L = L(:,order);

% Remove isolated nodes from the Facebook network.
bad = sum(A, 2) == 0;
A(bad, :) = [];
A(:, bad) = [];
L(bad, :) = [];
matches(bad) = [];

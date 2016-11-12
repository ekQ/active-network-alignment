function [A, B, L, matches] = get_genealogy_problem(problem_path, n, m)
% GET_GENEALOGY_PROBLEM Read a genealogy problem file.
%
% A: adjacency matrix of the source graph.
% B: adjacency matrix of the target graph.
% L: candidate matches and their unary similarities.
% matches: groundtruth
%
% Author: Eric Malmi (eric.malmi@gmail.com)

edges_A = load(strcat(problem_path, '/edges_A.txt'));
edges_A = edges_A + 1;  % Indexing starts from 0 in the original files.
A = sparse(n, n);
A(sub2ind(size(A), [edges_A(:,1); edges_A(:,2)], ...
          [edges_A(:,2); edges_A(:,1)])) = 1;


edges_B = load(strcat(problem_path, '/edges_B.txt'));
edges_B = edges_B + 1;  % Indexing starts from 0 in the original files.
B = sparse(m, m);
B(sub2ind(size(B), [edges_B(:,1); edges_B(:,2)], ...
          [edges_B(:,2); edges_B(:,1)])) = 1;

edges_L = load(strcat(problem_path, '/candidates_L.txt'));
edges_L(:,1:2) = edges_L(:,1:2) + 1;
L = sparse(n, m);
L(sub2ind(size(L), edges_L(:,1), edges_L(:,2))) = edges_L(:,3);

matches = (1:n)';
groundtruth = load(strcat(problem_path, '/groundtruth.txt'));
groundtruth = groundtruth + 1;
matches(groundtruth(:,1)) = groundtruth(:,2);

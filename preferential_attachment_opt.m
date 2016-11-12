% The original code is taken from: http://strategic.mit.edu/docs/matlab_networks/preferential_attachment.m
% Optimized by Eric Malmi

% Routine implementing a simple preferential attachment (B-A) model for network growth
% The probability that a new vertex attaches to a given old vertex is proportional to the (total) vertex degree
% Vertices arrive one at a time
% INPUTs: n - final (desired) number of vertices, m - # edges to attach at every step
% OUTPUTs: edge list, [number of edges x 3]
% NOTE: Assume undirected simple graph
% Source: "The Structure and Function of Complex Networks", M.E.J. Newman;  "Emergence of Scaling in Random Networks" B-A.
% GB, March 18, 2006

function E = preferential_attachment_opt(n,m)

E = zeros(2*m*n, 2);
A = sparse(n, n);
% start with an edge
E(1, :) = [1, 2];
E(2, :) = [2, 1];
Eidx = 3;
A(1,2) = 1;
A(2,1) = 1;

for vidx = 3:n
    deg = sum(A, 2);
    deg = full(deg(1:vidx-1,:));
    % add m edges
    r = randsampleWRW(1:(vidx-1), m, deg / sum(deg));
    for node=1:length(r)
      E(Eidx,:) = [r(node) vidx];
      Eidx = Eidx + 1;
      E(Eidx,:) = [vidx r(node)];
      Eidx = Eidx + 1;
      A(r(node), vidx) = 1;
      A(vidx, r(node)) = 1;
    end
end
% Remove extra rows.
E(Eidx:end,:) = [];
end

function v=randsampleWRW(x,k,w)
% same than randsample but Without Replacement and with Weighting

% Returns V, a weigthed sample of K elements taken among X without replacement
% X a vector of numerics
% K amount of element to sample from x
% W a vector of positive weights w, whose length is length(x)

% EXAMPLE:
% for i=1:100
% v(i)=randsampleWRW([0,0.5,3,20],1,[0.5,0.4,0.05,0.05]);
% end
% plot(v,'o')

% Was initialy made by ROY, i just correct some details and add help
% see the link:
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/141124
% Jean-Luc Dellis, april 2010

%if k>length(x), error('k must be smaller than length(x)'), end
%if ~isequal(length(x),length(w)),error('the weight W must have the length of X'),end
v=zeros(1,k);
for i=1:k
    v(i)=randsample(x,1,true,w);
    w(x==v(i))=0;
end
end

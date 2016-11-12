function [xbest,status,hist,mabest,mbbest,best_prob] = netalignmr(S,w,a,b,li,lj,gamma,stepm,rtype,maxiter,verbose)
% NETALIGNMR Compute the matching relaxation heuristic for network alignment
%
% Given a network alignment problem, the matching heuristic solves a
% sequence of matching problems to generate good upper and lower bounds on
% the solutions.
%
% [xbest,status,hist] = netalignmr(S,w,a,b,li,lj,stepm,rtype,
%    maxiter,verbose) 
% fully specifies all inputs and outputs.  
%
% Input
% -----
%   S : the complete set of squares
%   w : the matching weights for all edges in the link graph L
%   a : the value of alpha in the netalign objective
%   b : the value of beta in the netalign objective
%   li : the start point of each edge in L (a vertex number from graph A)
%   lj : the end point of each edge in L (a vertex number from graph B)
%   gamma : the starting step value (default = 0.5)
%   stepm : number of non-decreasing iterations adjusting the step length
%     (default = 100)
%   rtype : the rounding type (default = 1)
%     rtype = 1 : only consider the current matching
%     rtype = 2 : try enriching the matching with info from other squares
%   maxiter : maximum number of iterations to take (default = 1000)
%   verbose : output verbose information at each iteration (default = true)
%
% Output
% ------
%   xbest : the best heuristic solution, it may or may not be a matching
%   status : three values to describe the status of the solution
%     status(1) = 1 if the problem is solved to optimality, 0 otherwise
%     status(2) = best lower bound
%     status(3) = best upper bound
%   hist : the history of properties of each iteration
%     hist(k,) = 
%     hist(k,) = best lower bound at iteration k
%     hist(k,) = best upper bound at iteration k
%     hist(k,) = current value at iteration k
%     hist(k,) = matching weight at iteration k
%     hist(k,) = matching cardinality at iteration k
%     hist(k,) = overlap at iteration k
%
% Example:
%   load('../data/natalie_graphs');
%   netalignmr(S,w,0,1,li,lj);

% David F. Gleich, Ying Wang, and Mohsen Bayati
% Copyright, Stanford University, 2008-2009
% Computational Approaches to Digital Stewardship

% 2009-06-11: Initial coding (David and Ying)
% 2009-06-15: Cleanup and optimization (David)
% 2016-10-18: Return the best solution and the corresponding problem (Eric Malmi).

if ~exist('a','var') || isempty(a), a=1; end
if ~exist('b','var') || isempty(b), b=1; end
if ~exist('stepm','var') || isempty(stepm), stepm=25; end
if ~exist('rtype','var') || isempty(rtype), rtype=1; end
if ~exist('maxiter', 'var') || isempty(maxiter), maxiter=1000; end
if ~exist('verbose', 'var') || isempty(verbose), verbose=1; end 
if ~exist('gamma', 'var') || isempty(gamma), gamma = 0.4; end

m = max(li);
n = max(lj);

[rp ci ai tripi matn matm] = bipartite_matching_setup(w,li,lj,m,n); 
mperm = tripi(tripi>0); % a permutation for the matching problem

S = double(S);
U = sparse(size(S,1),size(S,2));

xbest = w;
xbest(:) = 0;

flower = 0;                       % best lower bound on the solution
fupper = Inf;                     % best upper bound on the solution
next_reduction_iteration = stepm; % reduce the step


hist = zeros(maxiter,7);

if verbose % print the header
    fprintf('%5s   %4s   %8s   %7s %7s %7s  %7s %7s %7s %7s\n', ...
        'best', 'iter', 'norm-u', 'lower','upper', 'cur', 'obj', 'weight', 'card', 'overlap');
end    
for iter = 1:maxiter   
    [q,SM] = maxrowmatch((b/2)*S + U-U',li,lj,m,n);
    x = a*w + q;
    
    %[val ma mb mi]= bipartite_matching(nw,li,lj,m,n);
    ai=zeros(length(tripi),1); ai(tripi>0)=x(mperm);
    [val ma mb mi]= bipartite_matching_primal_dual(rp,ci,ai,tripi,matn,matm);
    if verbose
        fprintf('Matching: n=%d, m=%d, nnz=%d, nnz_prop=%f\n', matn, matm, length(ci)-1, length(ci)/(matn*matm));
    end
        
    % compute statistics 
    matchval = mi'*w;
    overlap = mi'*S*mi/2; 
    card = length(ma);
    f = a*matchval + b*overlap;
    
    if val<fupper
        fupper=val;
        next_reduction_iteration = iter+stepm;
    end
    if f>flower
        flower=f;
        itermark = '*';
        xbest = mi;
        mabest = ma;
        mbbest = mb;
        problem = struct('rp', rp, 'ci', ci, 'ai', ai, 'tripi', tripi, ...
                'n', n, 'm', m);
        best_prob = problem;
    else
        itermark = ' ';
    end        
    
    if rtype==1
        % no work
    elseif rtype==2
        mw = S*x;
        mw = a*w + b/2*mw;
       
        ai=zeros(length(tripi),1); ai(tripi>0)=mw(mperm);
        [val ma mb mx] = bipartite_matching_primal_dual(rp,ci,ai,tripi,matn,matm);
        
        card = length(ma);
        matchval = mx'*w;
        overlap = mx'*S*mx/2;
        f = a*matchval + b*overlap;
        
        if f>flower
            flower=f;
            itermark = '**';
            mi = mx;
            xbest = mw;
        end
    end
    
    % report on current iter
    hist(iter,1:end) = [norm(nonzeros(U),1), flower, fupper, f, matchval, card, overlap];
    
    if verbose
        fprintf('%5s   %4i   %8.1e   %7g %7g %7g  %7g %7g %7i %7i\n', ...
            itermark, iter, norm(nonzeros(U),1), ...
            flower, fupper, val, ...
            f, matchval, card, overlap);
    end
    
    if iter==next_reduction_iteration
        gamma = gamma*0.5;
        if verbose
            fprintf('%5s   %4s   reducing step to %g\n', '', '', gamma);
        end
        if gamma < 1e-24, break; end
        next_reduction_iteration = iter+stepm;
    end
    
    if (fupper-flower)<1e-2
        break;
    end
    
    U = U - diag(sparse(gamma*mi))*triu(SM) + tril(SM)'*diag(sparse(gamma*mi));
    
    U = bound(U, -.5, .5);
end
hist = hist(1:iter,:);
status = zeros(1,3);
status(1) = (fupper-flower)<1e-2;
status(2) = flower;
status(3) = fupper;

function S=bound(S,a,b)
S=spfun(@(x) min(max(x,a),b), S);

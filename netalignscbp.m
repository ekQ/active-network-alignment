function [mbest hista histb Sk mabest mbbest best_prob] = netalignscbp(...
        S,w,a,b,li,lj,gamma,dtype,maxiter,verbose)
% NETALIGNMPSC Solve the network alignment problem with message passing
%
% This version of network alignment uses the matrix formulation of the 
% algorithm.

% David F. Gleich, Ying Wang, and Mohsen Bayati
% Copyright, Stanford University, 2010
% Computational Approaches to Digital Stewardship
%
% 2016-10-18: Return the best solution and the corresponding problem (Eric Malmi).

if ~exist('a','var') || isempty(a), a=1; end
if ~exist('b','var') || isempty(b), b=1; end
if ~exist('gamma','var') || isempty(gamma), gamma=0.99; end
if ~exist('dtype', 'var') || isempty(dtype), dtype=2; end
if ~exist('maxiter', 'var') || isempty(maxiter), maxiter=100; end
if ~exist('verbose', 'var') || isempty(verbose), verbose=1; end 

nedges = length(li); nsquares = nnz(S); m = max(li); n = max(lj);

% Initialize the messages
y = zeros(nedges,1); 
z = y; % copy z as the same size
Sk = zeros(nsquares, 1);  % messages and their "reversed/transposed" copies
Skt = zeros(nsquares, 1);
Sw = b * ones(nsquares, 1);
MSc = b * ones(nsquares, 4); MScnew = zeros(nsquares, 4); MP = zeros(nsquares, 4);
if dtype>1, d = y; end % needed for damping scheme
% Initialize a few parameters
damping = gamma; curdamp = 1; iter = 1; 

% Initialize history
hista = zeros(maxiter,4); % history of messages from ei->a vertices
histb = zeros(maxiter,4); % history of messages from ei->b vertices
fbest = 0; fbestiter = 0;
if verbose % print the header
fprintf('%4s   %4s   %7s %7s %7s %7s   %7s %7s %7s %7s\n', ...
        'best', 'iter', 'obj_ma', 'wght_ma', 'card_ma', 'over_ma', ...
        'obj_mb', 'wght_mb', 'card_mb', 'over_mb');
end

% setup the matching problem once
[rp ci ai tripi matn matm] = bipartite_matching_setup(...
						      w,li,lj,m,n);
mperm = tripi(tripi>0);
[rpS ciS] = sparse_to_csr(S);
riS = zeros(nsquares, 1);
for i=1:nedges
for j=rpS(i):rpS(i+1)-1
riS(j) = i;
    end
end

Stmp = sparse(riS, ciS, [1:nsquares]');
[ti tj transmap] = find(Stmp);

SC = findSC(S, li, lj);
r = 1;

while iter<=maxiter
    curdamp = damping*curdamp;
    if dtype>1, dold=d; end

    x = Sw;
    % updates for M_{ii'->f_i} and M_{ii'->g_{i'}}, y, z
    d = zeros(nedges, 1);
    Sk_tran = Sk(transmap);
    for i=1:nedges
        for j=rpS(i):rpS(i+1)-1
    %            d(i) = d(i) + max(0, b + Sk_tran(j)) - max(0, Sk_tran(j));
d(i) = d(i) + max(0, (r + 1) * Sw(j) + Sk_tran(j)) - max(0, r * Sw(j) + Sk_tran(j));
        end
    end
ynew = a*w - max(0,implicit_maxprod(m,li,z)) + d;
znew = a*w - max(0,implicit_maxprod(n,lj,y)) + d;

    % updates for M_{ii'->h_{ii'jj'}}, Sk
    t = 0;
    for i=1:nedges
        t = ynew(i) + znew(i) - a*w(i) - d(i);
        for j=rpS(i):rpS(i+1)-1
%            Skt(j) = t - max(0, b + Sk_tran(j)) + max(0, Sk_tran(j));
            Skt(j) = t - max(0, (r + 1) * Sw(j) + Sk_tran(j)) - max(0, r * Sw(j) + Sk_tran(j));
        end
    end
    
    St = Sk + Sk(transmap) + Sw;
    
    %updates for M_{ii'jj'->h_{ii'jj'}}, Sw
    for i=1:4
        MP(:,i) = max(0, implicit_maxprod(max(SC(:,i)), SC(:,i), MSc(:, i)));
    end
    Sw = -(MP(:,1) + MP(:,2) + MP(:,3) + MP(:,4)) + b;
    
    %updates for M_{ii'jj'->d_{ii'j}}, MSc
    tmp = Sk + Sk(transmap);
    for i=1:4
        MSc(:, i) = Sw + MP(:, i) + min(tmp, min(Sk, Sk(transmap)));
    end
    
    if dtype==1
        Sk = curdamp*Skt+ (1-curdamp)*Sk;
        y = curdamp*ynew+(1-curdamp)*y;
        z = curdamp*znew+(1-curdamp)*z;
    elseif dtype==2
        prev = (y+z-a*w+dold);
        y = ynew+(1-curdamp)*prev;
        z = znew+(1-curdamp)*prev;
        clear prev;
        Sk = Skt + (1-curdamp)*St;
        Sw = Sw + (1-curdamp)*St;
        for i=1:4
            MSc(:, i) = MSc(:, i) + (1-curdamp)*St;
        end
    elseif dtype==3
        prev = (y+z-a*w+dold);
        y = curdamp*ynew+(1-curdamp)*prev;
        z = curdamp*znew+(1-curdamp)*prev;
        clear prev;
        Sk = curdamp * Skt + (1-curdamp)*St;
        Sw = curdamp * Sw + (1-curdamp)*St;
        for i=1:4
            MSc(:, i) = curdamp * MSc(:, i) + (1-curdamp)*St;
        end
    end
    %tot = roundbyS(Sk, li, lj, riS, ciS, m, n);
    %tot2 = roundbyyz(S, y, z, li, lj, m, n);
    iter=iter+1;
    [hista(iter,:),ma,mb,problem] = round_messages(y,S,w,a,b,rp,ci,tripi,matn,matm,mperm);
    [histb(iter,:),ma2,mb2,problem2] = round_messages(z,S,w,a,b,rp,ci,tripi,matn,matm,mperm);

    if hista(iter,1)>fbest
        fbestiter=iter; mbest=y; fbest=hista(iter,1);
        mabest = ma;
        mbbest = mb;
        best_prob = problem;
    end
    if histb(iter,1)>fbest
        fbestiter=-iter; mbest=z; fbest=histb(iter,1);
        mabest = ma2;
        mbbest = mb2;
        best_prob = problem2;
    end
    %fbest = max(fbest, tot);
    %fbest = max(fbest, tot2);
    
    if verbose
        if fbestiter==iter, bestchar='*a'; 
        elseif fbestiter==-iter, bestchar='*b';
        else bestchar='';
        end
        fprintf('%4s   %4i   %7g %7g %7i %7i  %7g %7g %7i %7i\n', ...
            bestchar, iter, hista(iter,:), histb(iter,:));
    end
end

end

function S=bound(S,a,b)
  S=spfun(@(x) max(b+x,a)-max(x,a),S);
end

function y=maxprod(ai,aj,av,m,x) 
  y=-inf*ones(m,1);
  for i=1:numel(ai), y(ai(i)) = max(y(ai(i)),av(i)*x(aj(i))); end
end
function y=implicit_maxprod(n,ai,x)
% implicitly compute 
%   Ai=sparse(ai,1:length(ai),1,n,length(z))
     %   A=Ar'*Ar - speye(length(z))
%   maxprod(A,x)
% which has a very nice structure if you look at the actual summation
% definition.  It is just the maximum value 
  N = length(ai);
  y=-inf*ones(N,1);
  max1 = -inf*ones(n,1);
  max2 = -inf*ones(n,1);
  max1ind = zeros(n,1);
  for i=1:N
    if x(i)>max2(ai(i))
      if x(i)>max1(ai(i))
          max2(ai(i)) = max1(ai(i));
          max1(ai(i)) = x(i);
          max1ind(ai(i)) = i;
      else
          max2(ai(i)) = x(i);
      end
    end
  end
  for i=1:N
    if i==max1ind(ai(i))
      y(i)=max2(ai(i));
    else
      y(i)=max1(ai(i));
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [info,ma,mb,problem]=round_messages(messages,S,w,alpha,beta,rp,ci,tripi,n,m,perm)
ai=zeros(length(tripi),1); ai(tripi>0)=messages(perm);
[val ma mb mi]= bipartite_matching_primal_dual(rp,ci,ai,tripi,n,m);
problem = struct('rp', rp, 'ci', ci, 'ai', ai, 'tripi', tripi, 'n', n, 'm', m);
%fprintf('Matching: n=%d, m=%d, nnz=%d, nnz_prop=%f\n', n, m, length(ci)-1, length(ci)/(n*m));
matchweight = sum(w(mi)); cardinality = sum(mi); overlap = (mi'*(S*double(mi)))/2; 
f = alpha*matchweight + beta*overlap;
info = [f matchweight cardinality overlap];
end

function tot=roundbyS(Sk, li, lj, ri, ci, m, n)
    [v ind] = sort(Sk, 'descend');
    markA = zeros(m, 1);
    markB = zeros(n, 1);
    chosen = zeros(size(li, 1), 1);
    tot = 0;
    for i=1:size(Sk, 1)
        v1 = ri(ind(i));
        v2 = ci(ind(i));
        if (v1 > v2)
            continue;
        end
        if ((~chosen(v1) && (markA(li(v1)) || markB(lj(v1)))) || (~chosen(v2) && (markA(li(v2)) || markB(lj(v2)))))
            continue;
        end
        tot = tot + 1;
        chosen(v1) = 1;
        chosen(v2) = 1;
        markA(li(v1)) = 1;
        markA(li(v2)) = 1;
        markB(lj(v1)) = 1;
        markB(lj(v2)) = 1;
    end
end

function tot=roundbyyz(S, y, z, li, lj, m, n)
    [v ind] = sort(y+z, 'descend');
    markA = zeros(m, 1);
    markB = zeros(n, 1);
    x = zeros(size(y));
    for i=1:size(ind,1)
        if (~markA(li(ind(i))) && ~markB(lj(ind(i))))
            x(ind(i)) = 1;
            markA(li(ind(i))) = 1;
            markB(lj(ind(i))) = 1;
        end
    end
    tot = x' * S * x / 2;
end

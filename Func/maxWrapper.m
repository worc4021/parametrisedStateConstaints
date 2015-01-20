function sol = maxWrapper(P,q,Ch,phih,actW,actX,s,m)

% [K  L'][x]        =   [LHSconst]  + [LHSvar][theta]
% [L  0 ][lambda]   =   [LHSconst]  + [LHSvar][theta]

nlambda = length(find(actW{m}));
neta = length(find(actX{m}));

K = blkdiag(-s.gammasq*eye(s.nW),P(1:s.nX,1:s.nX));
L = [-s.D,eye(s.nX);    
    zeros(neta,s.nW),s.S{m}(actX{m},1:s.nX);
    s.G(actW{m},:),zeros(nlambda,s.nX);
    zeros(size(Ch,1),s.nW),Ch(:,1:s.nX)];
 
LHSconst = [zeros(s.nW,1);
            -q(1:s.nX);
            zeros(s.nX,1);
            s.s{m}(actX{m});
            ones(nlambda,1);
            phih];
            
LHSvar = [zeros(s.nW,s.nX+s.nU+s.nAlpha+s.nBeta);
          zeros(s.nX,s.nX+s.nU),-P(1:s.nX,s.nX+1:end);
          s.A,s.B,zeros(s.nX,s.nAlpha+s.nBeta);
          zeros(neta,s.nX+s.nU),-s.S(actX{m},s.nX+1:end);
          zeros(nlambda,s.nX+s.nU+s.nAlpha),ones(nlambda,s.nBeta);
          zeros(size(Ch,1),s.nX+s.nU),-Ch(:,s.nX+1:end)];
          
        [Q,~] = qr(L');
        if and(~isempty(L),size(L,1)+1<=size(Q,2))
            EW = eig(Q(:,size(L,1)+1:end)'*K*Q(:,size(L,1)+1:end));
        elseif size(L,1)+1>size(Q,2)
            EW = zeros(size(K));
        else
            EW = eig(K);
        end
        
        fprintf(' ')
        for i = 1:length(EW)
            fprintf('%6.1f ',EW(i))
        end
        fprintf(':')
      
[RHSconst,RHSvar,RHSorth] = KKTSolver(K,L,LHSconst,LHSvar);

firstIdx = 1;
Range = firstIdx:(firstIdx+s.nW-1);
firstIdx = Range(end)+1;

sol.wM = RHSvar(Range,:);
sol.wm = RHSconst(Range);

Range = firstIdx:(firstIdx+s.nX-1);
firstIdx = Range(end)+1;

sol.xpM = RHSvar(Range,:);
sol.xpm = RHSconst(Range);

Range = firstIdx:(firstIdx+s.nMu-1);
firstIdx = Range(end)+1;

sol.muM = RHSvar(Range,:);
sol.mum = RHSvar(Range);

if ~isempty(RHSorth)
    sol.muZ = RHSorth(Range,:);
else
    sol.muZ = [];
end

Range = zeros(1,0);

if neta > 0
    Range = firstIdx:(firstIdx+neta-1);
    firstIdx = Range(end)+1;
end

sol.etaM = RHSvar(Range,:);
sol.etam = RHSconst(Range);

if ~isempty(RHSorth)
    sol.etaZ = RHSorth(Range,:);
else
    sol.etaZ = [];
end

Range = zeros(1,0);

if nlambda > 0
    Range = firstIdx:(firstIdx+nlambda-1);
    firstIdx = Range(end)+1;
end

sol.lambdaM = RHSvar(Range,:);
sol.lambdam = RHSconst(Range);

if ~isempty(RHSorth)
    sol.lambdaZ = RHSorth(Range,:);
else
    sol.lambdaZ = [];
end

if ~isempty(Ch)
    sol.zetaM = RHSvar(firstIdx:end,:);
    sol.zetam = RHSconst(firstIdx:end);
    if ~isempty(RHSorth)
        sol.zetaZ = RHSorth(firstIdx:end,:);
    else
        sol.zetaZ = [];
    end
else
    sol.zetaM = [];
    sol.zetam = [];
    sol.zetaZ = [];
end

nHor = 2;

Lam = cell(1,nHor);
lam = cell(1,nHor);
xIdx = cell(1,nHor);
wIdx = cell(1,nHor);
uIdx = cell(1,nHor);

Lam{1} = Vout;
lam{1} = vout;


for  i = 1:nHor
    if i ~= 1
        [tLam,tlam] = stepBack(Lam{i-1},lam{i-1},sizes);
        [Lam{i},lam{i}] = PreConditioning(tLam,tlam);
    end
    xIdx{i} = false(1,size(Lam{i},1));
    wIdx{i} = false(1,size(G,1));
    uIdx{i} = false(1,size(F,1));
end

sizes.S = Lam;
sizes.s = lam;

figure(1)
hold on
for i = 1:nHor
    [Lt,lt] = rowReduce(Lam{i}(:,1:2),lam{i});
    plot(Polyhedron(Lt,lt),'color',[0,1,0],'alpha',.1)
end
hold off
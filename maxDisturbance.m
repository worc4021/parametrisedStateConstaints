% This script computes the maximal disturbance the system can cope with.

uMax = 3;

model.obj = [0,0,1];
model.A = sparse(Lam{1}(:,[1,2,4]));
model.sense = '<';
model.rhs = lam{1} - Lam{1}(:,3)*uMax;
model.modelsense = 'max';
model.varnames = {'x1', 'x2', 'beta'};

param.OutputFlag = 0;

res = gurobi(model,param);

clc;
fprintf('For alpha = %.3f the maximal beta is %.3f. \n',...
        uMax, res.objval);



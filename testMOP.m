FitnessFunction = @simple_multiobjective;
numberOfVariables = 1;
options = optimoptions(@gamultiobj,'PlotFcn',{@gaplotpareto,@gaplotscorediversity});
[x,fval] = gamultiobj(FitnessFunction,numberOfVariables,[],[],[],[],[],[],options);
function y = simple_multiobjective(x)
   y(1) = (x+2)^2 - 10;
   y(2) = (x-2)^2 + 20;
end


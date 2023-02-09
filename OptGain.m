
kp = 40;
kd = 2;
ki = 40;
alpha = 0.9;
b = 5.;
m1 = 10.;
m2 = 20.;
delta = 1.;


x0 = [0.21;0.21;1;1;1];
Optparameter = [kp,kd,ki,b,m1,m2,delta];
Optfun(Optparameter) 
FitnessFunction = @Optfun;
numberOfVariables = 7;
options = optimoptions("gamultiobj","MaxGenerations",10);
[x,fval] = gamultiobj(FitnessFunction,numberOfVariables,[],[],[],[],[],[],options);
function y = Optfun(x)
    kp = x(1);kd = x(2);ki = x(3);b = x(4);m1 = x(5);m2 = x(6);delta = x(7);
    Optparameter = [kp,kd,ki,b,m1,m2,delta];
    out = sim('pemfc_lumpekf2_prerun_2');
    
     y(1) = sum(abs(out.s.Data));
     y(2) = sum(abs(out.sdot.Data));
     y(3) = sum(abs(out.u.Data));
     y(4) = sum(abs(out.udot.Data));
     y(5) = sum(abs(out.u.Time.*out.s.Data.*out.s.Data));
end
 % Plot two objective functions on the same axis

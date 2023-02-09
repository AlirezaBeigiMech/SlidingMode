x0 = [0.21;0.21;0;0;0]; P0=diag([10,10,10000,10000,10000]);
Rsys=0.1^2;Qsys=diag([0.1^2 0.021^2 0.1^2/Area^2]);
xhat0 = [0;0;0;0;0]; 

rng(100)
dt=0.01;
t1=10;
tend=50;
%dt=0.05;tend=20;idx=1:5:Nsim; %Smaller
omega=2*pi*1;
if ( exist('J_ss','var') == 0 )
    J_ss = 10;
end
if ( exist('J_amp','var') == 0 )
    J_amp = J_ss/2;
end
if ( exist('xch_O2c','var') == 0 )
    xch_O2c = .21;
end
x0 = [xch_O2c;xch_O2c;J_ss;J_ss;J_ss]; 

Nsim=tend/dt;
t = (0:Nsim)*dt;
J_theo=J_ss+0*t;
ix_ = find(t>t1);J_theo(ix_) = J_ss+J_amp;

rng(100)
w = randn(Nsim+1,3)*Qsys.^0.5;
v = randn(Nsim+1,1)*Rsys.^0.5;
%xch_O2c=0.1;
if (add_noise==false)
    param.xch_H2a=xch_H2a+0*t;
    param.xch_O2c=xch_O2c+0*t+0*t;
    Jin=J_theo;v=v*0;
else
    param.xch_H2a=min(0.995, max(0,xch_H2a+w(:,1))');
    param.xch_O2c=max(0,xch_O2c+w(:,2)');
    Jin=max(0,J_theo+w(:,3)');
end
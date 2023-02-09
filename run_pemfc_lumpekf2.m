%function run_pemfc_lumpekf2(name,J_ss, xch_O2c, add_noise)

%%
%clear; 
close all
fno=1; f=0.5;

clearvars -except f;
name="sine"; J_ss = 1000; xch_O2c=0.21; add_noise=true;

sav=sprintf("output_%s_x%02dJ%04d_n%i",name,xch_O2c*100,J_ss,add_noise);
fprintf('Running %s',sav);
J_ss = 10;
name_ = name;load(sav);name=name_;
%t = (0:Nsim)*dt;

%%
param.xO2normal=0.21;
param.ekf.f=f;
param.sim.phiO2_eps = 1e-3;
param.sim.phieff_eps = 1e-3;

Qsys=diag([diag(Qsys);0.01*param.Ndot_a;0.01*param.Ndot_c]);
% Rsys=0.1^2;Qsys=diag([0.1^2 0.021^2 0.1^2/Area^2 0.01*param.Ndot_a 0.01*param.Ndot_c]);

%%
% nx=4;ny=1;ix=1;
% subplot(nx,ny,ix); plot(t,V+v); ix=ix+1;
% subplot(nx,ny,ix); plot(t,I); ix=ix+1;
% subplot(nx,ny,ix); plot(t,param.xch_H2a); ix=ix+1;
% subplot(nx,ny,ix); plot(t,param.xch_O2c); ix=ix+1;

simin.u.time = t';
%simin.u.signals.values=[param.xch_H2a' param.xch_O2c' Jin'];
disp('Changed u to \phi');
simin.u.signals.values=[param.Pa*xch_H2a+0*t' param.Pc*xch_O2c+0*t' J_theo' param.Ndot_a+0*t' param.Ndot_c+0*t'];
simin.u.signals.dimensions=5;

simin.y.time = t';
simin.y.signals.values=[V+v];
simin.y.signals.dimensions=1;

out = sim('pemfc_lumpekf2_prerun');
out.t= out.yhat(:,1);
idx1=find(27.2<out.t);idx2=find(out.t(idx1)<27.8);idx=idx1(idx2);
%%

set(gcf, 'PaperPositionMode', 'auto')   % Use screen size
plot(t,V,'k',out.yhat(:,1),out.yhat(:,2),'k--');
legend('System','EKF estimate');
xlabel('Time (s)'); ylabel('Voltage (V)');

%%
if param.ekf.f >= 0
    sav1=sprintf("_f%02d_%s_x%02dJ%04d_n%i",param.ekf.f*100,name,xch_O2c*100,J_ss,add_noise);
else
    sav1=sprintf("_f%s_%s_x%02dJ%04d_n%i","xx",name,xch_O2c*100,J_ss,add_noise);
end
sav=strcat("fig",sav1);
saveas(gcf, sav, 'png');saveas(gcf, sav, 'svg')

%%
V2ekf = interp1(t,V',out.yhat(:,1));

set(gcf, 'PaperPositionMode', 'auto')   % Use screen size
plot(out.yhat(:,1),out.yhat(:,2)-V2ekf,'k');
xlabel('Time (s)'); ylabel('Sys. Voltage - EKF estm. (V)');
%sav=sprintf("figerr2_%s_x%02dJ%04d_n%i",name,xch_O2c*100,J_ss,add_noise);

return

%%
SO2.Ecellc = SO.E0T+SO.Econc;

plot(t(idx),SO2.Ecellc([1 25 50],idx))
        
%% Compare plots
%out = [phi_eff,phi_avg,phiO2_out,Ecellc,Vohm,Vact,E0T,Econc];
n = size( out.mischat.signals.values );
misc_vals = reshape(out.mischat.signals.values,n(2),n(3));
phi_eff=misc_vals(1,:);phi_avg=misc_vals(2,:);phiO2_out=misc_vals(3,:);
Ecellc=misc_vals(4,:);Vohm=misc_vals(5,:);Vact=misc_vals(6,:);
E0T=misc_vals(7,:);Econc=misc_vals(8,:);
%Vcell = Ecellc - Vohm - Vact;
SO2.Ecellc = SO.E0T+SO.Econc;
fno=0;

%% Compare voltage
%figure(fno+1); plot(out.tout, out.yhat(:,2),t,V)
nsim=size(t,2); idx=1:10:nsim;
figure(fno+1); 
Nx=2;Ny=3;ix=1;
subplot(Ny,Nx,1); plot(out.tout, Ecellc, t(idx),SO2.Ecellc([1 25 50],idx));ylabel('Ecellc');
legend('EKF','01','25','50')
subplot(Ny,Nx,2); plot(out.tout, Vohm, t(idx),SO.Vohm([1 25 50],idx));ylabel('Vohm');
subplot(Ny,Nx,3); plot(out.tout, Vact, t(idx),SO.Vact([1 25 50],idx));ylabel('Vact');
subplot(Ny,Nx,4); plot(out.tout, E0T, t(idx),SO.E0T([1 25 50],idx));ylabel('E0T');
subplot(Ny,Nx,5); plot(out.tout, Econc, t(idx),SO.Econc([1 25 50],idx));ylabel('Econc');
subplot(Ny,Nx,6); plot(out.tout, out.yhat(:,2),t,V);ylabel('Econc');
%% Compare J/I -- sanity check ?
figure(fno+2); plot(out.tout,out.xhat(:,3), t,SO.JR([1 25 50],:))
figure(fno+3); plot(out.tout,out.xhat(:,3), t,SO.JR_lfilt([1 25 50],:))
figure(fno+4); plot(out.tout,out.xhat(:,4),t,SO.JR_lfilt1([1 25 50],:))

%% Compare phi effective
    %phi_eff = 0.1512*u(2) + 4/pi*x(1) - 4/(3*pi)*x(2) ...
    %    - (0.1894*u(3) + 8/pi^2*x(3))*0.21/param.J_lim;
figure(fno+5); 
Nx=2;Ny=4;ix=1;
subplot(Ny,Nx,ix);plot(out.tout,phi_eff, t(idx),SO.xO2c_conc([1 25 50],idx)); ix=ix+1; ylabel('\phi_{eff}')

subplot(Ny,Nx,ix);plot(simin.u.time,simin.u.signals.values(:,2), t(idx),xO2c([1 25 50],idx)); ix=ix+1;  ylabel('xO2c')
subplot(Ny,Nx,ix);plot(out.tout,out.xhat(:,1+1), t(idx),SO.xO2c_lfilt([1 25 50],idx)); ix=ix+1; ylabel('xO2c_{lfilt0}')
subplot(Ny,Nx,ix);plot(out.tout,out.xhat(:,1+2), t(idx),SO.xO2c_lfilt1([1 25 50],idx)); ix=ix+1; ylabel('xO2c_{lfilt1}')

subplot(Ny,Nx,ix);plot(simin.u.time(idx),simin.u.signals.values(idx,3), t(idx),J([1 25 50],idx)); ix=ix+1;  ylabel('J')
subplot(Ny,Nx,ix);plot(out.tout,out.xhat(:,1+5), t(idx),SO.JR([1 25 50],idx)); ix=ix+1;  ylabel('JR')
subplot(Ny,Nx,ix);plot(out.tout(idx),out.xhat(idx,1+3), t(idx),SO.JR_lfilt([1 25 50],idx)); ix=ix+1;ylabel('JR_{lfilt0}')




%%
%plot(out.tout, Vact,t,SO.Vact([1 25 50],:))
figure(fno+5); plot(out.tout, Vact,t,SO.Vact([1 25 50],:));
figure(fno+6); plot(out.tout, Vohm,t,SO.Vohm([1 25 50],:));

%%
idx1=find(27.6<out.t);idx2=find(out.t(idx1)<27.7);idx=idx1(idx2);
figure(fno+1); plot(out.t(idx),out.yhat(idx,2),'.')

%%
%    out = [phi_eff,phi_avg,phiO2_out,Ecellc,Vohm,Vact,E0T,Econc];

figure(fno+1);
n = size( out.mischat.signals.values );
misc_vals = reshape(out.mischat.signals.values,n(2),n(3));

n=n(1)*n(2); nr=ceil((n-1)/2);
for i=1:n
    subplot(nr,2,i); plot(out.t(idx),misc_vals(i,idx))
end

%%
p = plot(t,[t*0+xch_O2c; xO2c([25,50],:)]);
for i=1:3
    p(i).Color=[0 0 0];
end
p(2).LineStyle="--";p(3).LineStyle=":";
xlabel('Time(s)'); ylabel('O2 Mole fraction');ylim([0 0.35])
legend('Inlet','Middle','Outlet')

%%
p = plot(t,[t*0+xch_H2a; xH2a([25,50],:)]);
for i=1:3
    p(i).Color=[0 0 0];
end
p(2).LineStyle="--";p(3).LineStyle=":";
xlabel('Time(s)'); ylabel('H2 Mole fraction');ylim([0.95 1.05])
legend('Inlet','Middle','Outlet')
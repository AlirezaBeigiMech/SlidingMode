%%Default
%%
%        output[i,0] = np.sum(np.absolute(np.diff(BATCH_u))) 
%        output[i,1] = np.sum(np.absolute(BATCH_S1))
%        output[i,2] = np.sum(np.absolute(BATCH_S2))
%        output[i,3] = np.sum(np.absolute(np.dot(BATCH_t,np.dot(BATCH_S1,BATCH_S1))))
%        output[i,4] = np.sum(np.absolute((BATCH_u)))
%%	
la = 100e-2; lc = la; w = 5e-2; slm = 1/24.4/60; %slm to moles/sec;
    Area = la*w;

	%Universal Constantsin %SI Units
	param.R=8.314; param.F=96485; param.DH2OH2=0.915e-4; param.DH2OO2=0.282e-4; 
	param.DN2O2=0.220e-4; param.DCO2O2=0.150e-4;

	%FC Reaction parameters (Per cell values from Wang et al.)

	param.act_b = 1.0417e-04; param.J0 = 0.005;	
	param.E0_cell = 1.2271; param.kE=1.7708e-05;
	param.eta_0 = 0.4197 - 298*param.act_b*(la*w/0.34747);
	param.act_a = -0.0029  - param.act_b*(la*w/0.34747);
	%param.J_lim = 2.6e4; %Limiting current density


	param.la = la; param.lc = lc; param.w = w;

	param.tau_flux_O2 = 0.1;
	param.Cdbl = 4.8000; %Double layer Cap
	param.Rther = 2.7778;
	param.Cther = 0.7200;

	param.r1 = 0.006;
	param.r2 = 4.0000e-05;
	param.rt = 5.0000e-05;

	%----------Rel---
	param.xch_H2a=0.99;
	param.xch_H2c=0;
	xch_O2c=0.21;

	param.Pa = 1.5;	%In atm
	param.Pc = 1; %atm.
	param.Tfc = 273+33;

	%--------
	param.Ndot_a = 1.03*slm; %Moles per sec of non particapating gas
	param.Ndot_c = 1.03*slm; %Moles per sec of non particapating gas
	
	%param.r_wat = 0.1; param.r_wat_max = 0.1;
	
	%xO2c_conc = xO2c[i-1]*(1-J/J_lim) - 0.10*( 1-exp(-0.1*J) );
	param.r_dO2_prexp = 0;
	param.r_dO2_exp = 0.1;
	param.r_dO2_Jmin = 200;
	param.r_dO2_tau = param.tau_flux_O2;
	
	%param.cool_Tin = 273+68.5;
	%param.cool_Tout = 273+68.5;
	param.Tfc =  273+68.5;
	%param.Troom = 298;
	
	param.Bconc = 0.05;
	
	param.exp_ratio = 2;

    %param.J_lim = 304.24 * xch_O2c;
   
%%Parameters: Updated from email Mar 14 '14

%cathode (air)

%Dont change:
slm = 1/24.4/60; %slm to moles/sec
Nele = 50; 

%Change below:

%Known:
%Anode H2:
xch_H2a=1-0.01;
%Cathode changes:
%Set ambient and inital flow
xch_N2c_amb=0.78; xch_H2Oc_amb = 1/100; xch_CO2c_amb=0.038/100; 
xch_O2c_amb = 1 - (xch_N2c_amb + xch_CO2c_amb + xch_H2Oc_amb );

%Have some idea:
param.J_lim = 3e4;	%Not sure 2.6 - 3.6 A/cm^2
param.xO2normal=0.21;

%No idea
param.r_dO2_prexp = 0; param.r_dO2_exp  = 0.004;
%param.r_dO2_prexp = 0.0375; param.r_dO2_exp  = 0.008;

%Need to set this in others:
%param.act_b = 1.0417e-04;
%param.eta_0 = 0.4197 - 298*param.act_b*(param.la*param.w/0.34747);
%param.eta_0 = 0.2197 - 298*param.act_b*(param.la*param.w/0.34747);

param.tau_flux_O2 = 6;

flag_start = 8;
flag_dyn = 4;

%%Calculations
eta_0 = param.eta_0 + param.act_a*(param.Tfc-298);
param.BV_Vact0 = param.act_b*param.Tfc/2; param.BV_I0 = exp( -eta_0/(2*param.BV_Vact0) );
param.BV_J0 = param.BV_I0/(param.la*param.w);
param.BV_VactH20 = 8.615e-5; param.BV_JH20 = 1e2;

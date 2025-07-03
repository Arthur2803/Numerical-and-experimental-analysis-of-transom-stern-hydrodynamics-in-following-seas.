
clc
clear all
close all

% Ship velocities in m/s
Vsvec = [0.1:0.01:1.3];
L = 1.06;                    % Model length
Lwl = 1;                  % Waterline length
Los = 1.03;                  % Overall ship length
T = 0.01;                  % Draft
B = 0.25;                   % Beam
S = 0.3159;                 % Wetted surface area
CB = 0.703;
TA = T;
TF = T;
Dp = 0.001;
NRud = 0;
NBrac = 0;
NThr = 0;
NBoss = 0;
k = 0.3175;		% Form factor to the ship

rho = 1000;
gravk = 9.81;   %Gravity
nu = 1.1395E-6; % viscosity


%Calculation of 'Froude length', Lfn:
if Los/L < 1
   Lfn = Los;
elseif (Los/L >= 1) && (Los/L < 1.1)
   Lfn = L+2/3*(Los-L);
elseif Los/L >= 1.1
   Lfn = 1.0667*L;   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants from Hollenbachs paper: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 'Mean' resistance coefficients
   a = [-0.3382 0.8086 -6.0258 -3.5632 9.4405 0.0146 0 0 0 0];  	
%a1 means a(1) and so on
   b = [-0.57424 	 13.3893	90.5960; 	%b12 means b(1,2)
	   4.6614	-39.721	-351.483;
      -   1.14215	-12.3296	459.254];
   d = [0.854 -1.228 0.497];
   e = [2.1701 -0.1602];
   f = [0.17 0.20 0.60];
   g = [0.642 -0.635 0.150];

   % 'Minimum' resistance coefficients
   a_min = [-0.3382 0.8086 -6.0258 -3.5632 0 0 0 0 0 0];
   b_min = [-0.91424 13.3893 90.5960;...
         4.6614 -39.721 -351.483;...
         -1.14215 -12.3296 459.254];
   d_min = [0 0 0];
   e_min = [1 0];
   f_min = [0.17 0.2 0.6];
   g_min = [0.614 -0.717 0.261];


cc = 0;
% Loop over velocities
for Vs = Vsvec
	
	cc = cc + 1;
	
	% Froude's number
	Fn = Vs/sqrt(gravk*Lfn);			
	Fnkrit = d*[1 CB CB^2]';
    
	c1 = Fn/Fnkrit;
	c1_min = Fn/Fnkrit;
	
	Rns = Vs*L/nu;						% Reynold's number for ship
	CFs = 0.075/(log10(Rns)-2)^2;			% ITTC friction line for ship				

	% Calculation of C_R for given ship 
	% Mean value
	
	CRFnkrit = max(1.0,(Fn/Fnkrit)^c1);
	
	kL = e(1)*L^(e(2));
	
	% There is an error in the hollenbach paper and in Minsaas' 2003 textbook, which
	% is corrected in this formula by dividing by 10
	CRstandard = [1 CB CB^2]*(b*[1 Fn Fn^2]')/10;
	
	CR_hollenbach = CRstandard*CRFnkrit*kL*prod([T/B B/L Los/Lwl Lwl/L (1+(TA-TF)/L) ...
			Dp/TA (1+NRud) (1+NBrac) (1+NBoss) (1+NThr)].^a);
	
	CR = CR_hollenbach*B*T/S;   			% Resistance coefficient, scaled for wetted surface
	C_Ts = CFs + CR;				% Total resistance coeff. ship 
	R_T_mean = C_Ts*rho/2*Vs^2*S;			% Total resistance to the ship

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % When accounting for form factor, roughness and correlation coef.,  
    % as given by Minsaas, but using model‐scale friction:
    
    Rnm    = Vs * L / nu;                                 % model Reynolds no.
    CFm    = 0.075/(log10(Rnm)-2)^2;                      % model ITTC friction
    raw = 110.31*(150*Vs/0.514)^0.21 - 403.33;
    dCF = max(0, raw) * CFm^2;
    % roughness based on CFm
    CA     = -0.228e-3;                                   % correlation coefficient
    
    CR_2   = CR_hollenbach * B * T / S - k * CFm;         % corrected residual
    C_Ts_2 = (1 + k) * (CFm + dCF) + CR_2 + CA;           % model‐scale total C_T
    
    R_T_mean_2 = C_Ts_2 * rho/2 * Vs^2 * S;               % Total resistance [N]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	% Minimum values

	% There is an error in the hollenbach paper and in Minsaas' 2003 textbook, which
	% is corrected in this formula by dividing by 10
	CRstandard_min = [1 CB CB^2]*(b_min*[1 Fn Fn^2]')/10;
	
	CR_hollenbach_min = CRstandard_min*prod([T/B B/L Los/Lwl Lwl/L (1+(TA-TF)/L) ...
			Dp/TA (1+NRud) (1+NBrac) (1+NBoss) (1+NThr)].^a_min);
	
	CR_min = CR_hollenbach_min*B*T/S;
	
	% Total resistance coefficient of the ship 
	C_Ts_min = CFs + CR_min;				
	% Total resistance	
	R_T_min = C_Ts_min*rho/2*Vs^2*S;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% When accounting for form factor, roughness and correlation coef.,  as given by Minsaas
	CR_min_2 = CR_hollenbach_min*B*T/S - k*CFm;		% Resistance coefficient 
	C_Ts_min_2 = (1+k)*(CFs + dCF) + CR_min_2 + CA;	% Total resistance coeff. ship 
	R_T_min_2 = C_Ts_min_2*rho/2*Vs^2*S;		% Total resistance to the ship
		
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
    % Efficiency factors
    eta_D = 0.70; % Propulsion efficiency
    eta_M = 0.97; % Mechanical efficiency
	eta_tot = eta_D*eta_M;
	
	% Delivered power on shaft
	P_B_mean_2 = R_T_mean_2*Vs/eta_tot/1000; 		%[kW]
	P_B_min_2 = R_T_min_2*Vs/eta_tot/1000; 		%[kW]

	% Store results for plotting
	CFsvec(cc) = CFs;
	CRvec(cc) = CR;
	C_Tsvec(cc) = C_Ts;
	C_Ts_2vec(cc) = C_Ts_2;
	R_T_meanvec(cc) = R_T_mean;
	R_T_mean_2vec(cc) = R_T_mean_2;
	CR_minvec(cc) = CR_min;
	C_Ts_minvec(cc) = C_Ts_min;
	C_Ts_min_2vec(cc) = C_Ts_min_2;
	R_T_minvec(cc) = R_T_min;
	R_T_min_2vec(cc) = R_T_min_2;
	P_B_mean_2vec(cc) = P_B_mean_2;
	P_B_min_2vec(cc) = P_B_min_2;
		
end


% or, as CSV if you prefer:
T = table(Vsvec', C_Tsvec', C_Ts_2vec', ...
    'VariableNames',{'Vs','CTs','CTs2'});
writetable(T,'hollenbach_resultsNy.csv');


figure
plot(Vsvec,CFsvec,'k')
hold on
plot(Vsvec,CRvec,'k:')
plot(Vsvec,C_Tsvec,'k-.')
plot(Vsvec,C_Ts_2vec,'k--')
legend('CFs','CR','CTs','CTs2')
title('Mean values')
xlabel('Speed [m/s]')

figure
plot(Vsvec,CFsvec,'k')
hold on
plot(Vsvec,CR_minvec,'k:')
plot(Vsvec,C_Ts_minvec,'k-.')
plot(Vsvec,C_Ts_min_2vec,'k--')
legend('CFs','CR','CTs','CTs2')
title('Min values')	
xlabel('Speed [m/s]')

figure
plot(Vsvec,P_B_mean_2vec,'k')
hold on
plot(Vsvec,P_B_min_2vec,'k-.')
legend('Mean','Min')
title('Delivered power [kW]')	
xlabel('Speed [m/s]')

disp(C_Tsvec)
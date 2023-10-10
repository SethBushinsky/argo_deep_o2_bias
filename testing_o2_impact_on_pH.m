% investigating o2 impact on pH:

o2_offset = -20:10:20;
o2_mag = 300; %50:50:300;
no3_offset = -0.25:.25:0.25;
no3_mag = 5;
pH_out = NaN(length(o2_offset), length(o2_mag));
TA_out = NaN(length(o2_offset), length(o2_mag));

MeasIDVec = [1, 7, 3, 6];

for m=1:length(no3_offset)
    for o = 1:length(o2_offset)
        
        Coordinates_all = [-149.968     ,  -52.543     , 1500.        ];
        Measurements_all = [34.50930023,   2.80669999,      no3_mag - no3_offset(m),   o2_mag(1) - o2_offset(o)];
%         Measurements_offset = [ 34.50930023,   2.80669999,           172.22309662];

        pH_out(o,m) = LIPHR(Coordinates_all, Measurements_all, MeasIDVec, 'OAAdjustTF', 0, 'verboseTF', 0);
        TA_out(o,m) = LIAR(Coordinates_all, Measurements_all, MeasIDVec, 'verboseTF', 1);

    end
end

%% NO3 and Oxygen on pH at one point
MeasIDVec = [1, 7, 6];
MeasIDVec_ESPER = [1, 2, 6];

Coordinates_all = [-149.968     ,  -52.543     , 1500.        ];
Measurements_all = [34.50930023,   2.80669999,     300];

DesiredVariables = 3;

[Estimates_orig,Uncertainties]=ESPER_Mixed(DesiredVariables,Coordinates_all,Measurements_all,MeasIDVec_ESPER,'Equations', 7, 'EstDates', 2020, 'VerboseTF', 0);


%%



orig_pH = LIPHR(Coordinates_all, Measurements_all, MeasIDVec, 'OAAdjustTF', 0, 'verboseTF', 0);
% orig_TALK = LIAR(Coordinates_all, Measurements_all, MeasIDVec, 'verboseTF', 0);
disp('PH, Alk')
disp(['LPIHR ' num2str(orig_pH)])
% disp(orig_TALK)
Measurements_all = [34.50930023,   2.80669999,      300-5];

disp('Change oxygen by -5')
pH_O2 = LIPHR(Coordinates_all, Measurements_all, MeasIDVec, 'OAAdjustTF', 0, 'verboseTF', 0);
[Estimates_o2_change,Uncertainties]=ESPER_Mixed(DesiredVariables,Coordinates_all,Measurements_all,MeasIDVec_ESPER,'Equations', [7 15 16], 'EstDates', 2020);

disp('pH change LIPHR')
pH_O2 - orig_pH

% TA_O2 = LIAR(Coordinates_all, Measurements_all, MeasIDVec, 'verboseTF', 0);
% 
% disp('TA change')
% TA_O2 - orig_TALK
% 
% Measurements_all = [34.50930023,   2.80669999,       300];
% 
% pH_no3 = LIPHR(Coordinates_all, Measurements_all, MeasIDVec, 'OAAdjustTF', 0, 'verboseTF', 0);
% pH_no3 - orig_pH
% 
% Measurements_all = [34.50930023,   2.80669999,     300-5];
% 
% pH_O2_no3 = LIPHR(Coordinates_all, Measurements_all, MeasIDVec, 'OAAdjustTF', 0, 'verboseTF', 0);
% pH_O2_no3 - orig_pH

%% 

clf
subplot(1,1,1); hold on
for m = 1:length(o2_mag)

    plot(o2_offset, pH_out(:,m) - pH_out(3,m))

end

%% testing pH change of respiration on pCO2 - take surface water and progressively change DIC and Alk

% start with surface ocean values for DIC and Alk - say in SAMW - values
% from Carter et al. 2014
model_time = 250;

TA = NaN(model_time,1);
TA(1) = 2273.1;
DIC = NaN(model_time,1);
DIC(1) = 2122.2;
SAL = 34.181;
Theta = 5.247;
i=1;
Press_in = 0;
SI = 4.878;
PO4 = 1.44;

pH = NaN(model_time,1);
pCO2 = NaN(model_time,1);
O2 = NaN(model_time,1);
Revelle = NaN(model_time,1);
pH_Revelle = NaN(model_time,1);

% 1 mol of oxygen change decrease from respiration:
o2_change = -1;
delta_DIC = o2_change.*-106/154;
delta_ALK = delta_DIC*-16/106;
O2(1) = GGo2_units(Theta, SAL, 'umol');

for i = 1:model_time-1
    [DATA,~,~]= CO2SYSSOCCOM_smb(TA(i), DIC(i) , ...
        1,2,SAL, Theta, ...
        Theta ,...
        Press_in,Press_in,SI, PO4,1,10,3);

    TA(i+1) = TA(i) + delta_ALK;
    DIC(i+1) = DIC(i) + delta_DIC;
    O2(i+1) = O2(i) + o2_change;

    pH(i) = DATA(3);
    pCO2(i) = DATA(4);

    [DATA,~,~]= CO2SYSSOCCOM_smb(TA(i), DIC(i)+1, ...
        1,2,SAL, Theta, ...
        Theta ,...
        Press_in,Press_in,SI, PO4,1,10,3);

    pCO2_change_1_umol_DIC = DATA(4);
    pH_change_1_umol_DIC = DATA(3);
    Revelle(i) = pCO2_change_1_umol_DIC - pCO2(i);
    pH_Revelle(i) = pH_change_1_umol_DIC - pH(i);
end
%%
i=1
[DATA,~,~]= CO2SYSSOCCOM_smb(TA(i), DIC(i) , ...
        1,2,SAL, Theta, ...
        Theta ,...
        Press_in,Press_in,SI, PO4,1,10,3);
DATA_orig = DATA;

[DATA,~,~]= CO2SYSSOCCOM_smb(TA(i)-1, DIC(i) , ...
        1,2,SAL, Theta, ...
        Theta ,...
        Press_in,Press_in,SI, PO4,1,10,3);

DATA - DATA_orig
%%
[DATA_orig,~,~]= CO2SYSSOCCOM_smb(2337.788, 8.085734 , ...
        1,3,SAL, Theta, ...
        Theta ,...
        Press_in,Press_in,SI, PO4,1,10,3);

[DATA,~,~]= CO2SYSSOCCOM_smb(2337.788-.8, 8.085734 , ...
        1,3,SAL, Theta, ...
        Theta ,...
        Press_in,Press_in,SI, PO4,1,10,3);


DATA - DATA_orig

%%

plot_dir =     '/Users/smb-uh/UHM_Ocean_BGC_Group Dropbox/Seth Bushinsky/Work/Presentations/2023_04 BBB/';


clf
subplot(2,2,1)
plot(DIC, 'k')
ylabel('DIC')
xlabel('\mumol O2 respired')
yyaxis right
plot(TA)
ylabel('TA')

subplot(2,2,2)
plot(pH, 'k')
ylabel('pH')
yyaxis right
plot(diff(pH)*1000)
ylabel('\Delta mpH for 1 \mumol/kg O2 respiration')
xlabel('\mumol O2 respired')

subplot(2,2,3)
plot(O2, 'k')
ylabel('Oxygen')
xlabel('\mumol O2 respired')

subplot(2,2,4)
plot(Revelle, 'k')
ylabel('Revelle Factor')
yyaxis right
plot(pH_Revelle*1000)
ylabel('\Delta mpH / \Delta DIC')
xlabel('\mumol O2 respired')

print(gcf, '-dpng', '-r800',  [plot_dir  'Respiration_experiment_v2.png'])

%% plots checking pCO2 calculations in derived_3 files

clf; subplot(1,3,1); pcolor(TALK_LIAR); shading flat; colorbar; subplot(1,3,2); pcolor(TALK_w_NO3_ADJUST); ; shading flat; colorbar; 
subplot(1,3,3); pcolor(TALK_w_NO3_ADJUST - TALK_LIAR); shading flat; colorbar


%%
set(gcf, 'colormap', flipud(brewermap(20, 'Spectral')))
NO3_ADJUSTED = ncread('1902304_derived_3.nc', 'NITRATE_ADJUSTED');
NITRATE_ADJUSTED_w_O2_ADJUST = ncread('1902304_derived_3.nc', 'NITRATE_ADJUSTED_w_O2_ADJUST');

TALK_LIAR = ncread('1902304_derived_3.nc', 'TALK_LIAR');
TALK_LIAR_recalc = ncread('1902304_derived_3.nc', 'TALK_LIAR_recalc');

TALK_w_NO3_ADJUST = ncread('1902304_derived_3.nc', 'TALK_w_NO3_ADJUST');

DOXY_ADJUSTED = ncread('1902304_derived_3.nc', 'DOXY_ADJUSTED');
%%
clf; 

subplot(4,3,1); pcolor(juld_grid, press, DOXY_ADJUSTED); shading flat; colorbar; set(gca, 'ydir', 'reverse')
subplot(4,3,2); pcolor(juld_grid, press, DOXY_ADJUSTED);  shading flat; colorbar; set(gca, 'ydir', 'reverse')
subplot(4,3,3); pcolor(juld_grid, press, DOXY_ADJUSTED - DOXY_ADJUSTED); shading flat; colorbar; set(gca, 'ydir', 'reverse')

subplot(4,3,4); pcolor(juld_grid, press, NO3_ADJUSTED); shading flat; colorbar; set(gca, 'ydir', 'reverse')
subplot(4,3,5); pcolor(juld_grid, press, NITRATE_ADJUSTED_w_O2_ADJUST);  shading flat; colorbar; set(gca, 'ydir', 'reverse')
subplot(4,3,6); pcolor(juld_grid, press, NITRATE_ADJUSTED_w_O2_ADJUST - NO3_ADJUSTED); shading flat; colorbar; set(gca, 'ydir', 'reverse')

subplot(4,3,7); pcolor(juld_grid, press, TALK_LIAR); shading flat; colorbar; set(gca, 'ydir', 'reverse')
subplot(4,3,8); pcolor(juld_grid, press, TALK_LIAR_recalc);  shading flat; colorbar; set(gca, 'ydir', 'reverse')
subplot(4,3,9); pcolor(juld_grid, press, TALK_LIAR_recalc - TALK_LIAR); shading flat; colorbar; set(gca, 'ydir', 'reverse'); caxis([-2 2])

subplot(4,3,10); pcolor(juld_grid, press, pCO2_pH_orig_TALK_orig); shading flat; colorbar; set(gca, 'ydir', 'reverse')
subplot(4,3,11); pcolor(juld_grid, press, pCO2_pH_orig_TALK_no3);  shading flat; colorbar; set(gca, 'ydir', 'reverse')
subplot(4,3,12); pcolor(juld_grid, press, pCO2_pH_orig_TALK_no3 - pCO2_pH_orig_TALK_orig); shading flat; colorbar; set(gca, 'ydir', 'reverse'); %caxis([-2 2])


% investigating o2 impact on pH:

o2_offset = -20:10:20;
o2_mag = 300; %50:50:300;
no3_offset = -0.25:.25:0.25;
no3_mag = 5;
pH_out = NaN(length(o2_offset), length(o2_mag));
MeasIDVec = [1, 7, 3, 6];

for m=1:length(o2_mag)
    for o = 1:length(o2_offset)
        Coordinates_all = [-149.968     ,  -52.543     , 1500.        ];
        Measurements_all = [34.50930023,   2.80669999,      no3_mag,   o2_mag(m) - o2_offset(o)];
%         Measurements_offset = [ 34.50930023,   2.80669999,           172.22309662];

        pH_out(o,m) = LIPHR(Coordinates_all, Measurements_all, MeasIDVec, 'OAAdjustTF', 0, 'verboseTF', 0);
    end
end

%% NO3 and Oxygen on pH at one point
MeasIDVec = [1, 7, 3, 6];

Coordinates_all = [-149.968     ,  -52.543     , 1500.        ];
Measurements_all = [34.50930023,   2.80669999,      4,   300];
orig_pH = LIPHR(Coordinates_all, Measurements_all, MeasIDVec, 'OAAdjustTF', 0, 'verboseTF', 0);

Measurements_all = [34.50930023,   2.80669999,      4,   300-5];

pH_O2 = LIPHR(Coordinates_all, Measurements_all, MeasIDVec, 'OAAdjustTF', 0, 'verboseTF', 0);
pH_O2 - orig_pH

Measurements_all = [34.50930023,   2.80669999,      4+5*16/154,   300];

pH_no3 = LIPHR(Coordinates_all, Measurements_all, MeasIDVec, 'OAAdjustTF', 0, 'verboseTF', 0);
pH_no3 - orig_pH

Measurements_all = [34.50930023,   2.80669999,      4+5*16/154,   300-5];

pH_O2_no3 = LIPHR(Coordinates_all, Measurements_all, MeasIDVec, 'OAAdjustTF', 0, 'verboseTF', 0);
pH_O2_no3 - orig_pH

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

[DATA,~,~]= CO2SYSSOCCOM_smb(TA(i), DIC(i) , ...
        1,2,SAL, Theta, ...
        Theta ,...
        Press_in,Press_in,SI, PO4,1,10,3);
DATA_orig = DATA;

[DATA,~,~]= CO2SYSSOCCOM_smb(TA(i)-.15, DIC(i)+1 , ...
        1,2,SAL, Theta, ...
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
plot(diff(pH)*100)
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

print(gcf, '-dpng', '-r800',  [plot_dir  'Respiration_experiment.png'])

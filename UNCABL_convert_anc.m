%% Analyze Forces from UNCABL
% Don't know who's original code

clear
clc
% 
%        files = {'Fig8_1.anc';'Fig8_2.anc';'Fig8_3.anc';'Fig8_4.anc';'Fig8_5.anc';'GaitInit_L_1.anc';'GaitInit_L_2.anc';'GaitInit_L_3.anc';'GaitInit_L_4.anc';'GaitInit_L_5.anc';...
%           'GaitInit_R_1.anc';'GaitInit_R_2.anc';'GaitInit_R_3.anc';'GaitInit_R_4.anc';'GaitInit_R_5.anc';'OGFast_1.anc';'OGFast_2.anc';'OGFast_3.anc';'OGFast_4.anc';'OGFast_5.anc';'OGFast_6.anc'...
%         ;'OGPref_1.anc';'OGPref_2.anc';'OGPref_3.anc';'OGPref_4.anc';'OGPref_5.anc';'OGPref_6.anc';'OGSlow_1.anc';'OGSlow_2.anc';'OGSlow_3.anc';'OGSlow_4.anc';'OGSlow_5.anc';'OGSlow_6.anc';...
%           'PrecisionStepping_1.anc'; 'TMFast_1.anc'; 'TMPref_1.anc';'TMSlow_1.anc'}; %add all trials here
   
   
files = {'LFP_Left_1.anc';'LFP_Right_1.anc';'Norm_1.anc'; 'Slip_1.anc'; 'VR20_1.anc'; 'VR35_1.anc'}; %add all trials here

%I just run one subject at a time mostly
subjects = [{'OANF10\Session2'}];

FPcal_file = 'D:\RESEARCH\Projects\Codes\Codes used in Per Motor Rep\forcepla.cal'; % Lab Computer

for l=1:length(subjects)
    subject = char(subjects(l))
    odir = strcat(['D:\RESEARCH\Projects\Data Collections\R21Repertoire\YA\clean\', subjects{l},'\Generated_Files\']); % Update to your subject's Cortex Files
 
for t = 1:length(files)
    input_file = char(files(t))
    [mass(t)] = convertFPdata(input_file,FPcal_file,odir);

end
% submass(l) = mean(mass);
end

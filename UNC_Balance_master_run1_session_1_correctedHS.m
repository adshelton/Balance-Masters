%Master analayis routine 1 for walking data for balance analyses
% Uses kinematics to define gait cycle events

% Input: 
%       Cortex "Generated files" (.trc and .anc)
%       Treadmill Speed
% Output: 
%       Step kinematics (length and width)
%       Time series of EMG data (if available)
%       Time series of select markers of interest (i.e., C7, S2, etc.)


clear all
close all
clc
addpath('D:\RESEARCH\Projects')

%ADS edited on 17 May 2021 for Per Motor Rep
%ADS edited on 23 June 2022 to remove old commented out code for SW and SL

%% Settings
sn='OANF09\Session1'; %Subject Number YA00\Session1
fs=100; %sampling frequency

% Define the path where your data is located:
path=char(['D:\RESEARCH\Projects\Data Collections\R21Repertoire\YA\clean\' sn '\Generated_files\']); %you will need to customized this.%E:\RESEARCH\Projects\Data Collections\R21Repertoire\YA\YA03\Session1

%Select Trials of Interest to Compile (list all as needed):

 files = {'Fig8_1';'Fig8_2';'Fig8_3';'Fig8_4';'Fig8_5';'GaitInit_L_1';'GaitInit_L_2';'GaitInit_L_3';'GaitInit_L_4';'GaitInit_L_5';...
    'GaitInit_R_1';'GaitInit_R_2';'GaitInit_R_3';'GaitInit_R_4';'GaitInit_R_5';'OGFast_1';'OGFast_2';'OGFast_3';'OGFast_4';'OGFast_5';'OGFast_6'...
   ;'OGPref_1';'OGPref_2';'OGPref_3';'OGPref_4';'OGPref_5';'OGPref_6';'OGSlow_1';'OGSlow_2';'OGSlow_3';'OGSlow_4';'OGSlow_5';'OGSlow_6';...
    'PrecisionStepping_1'; 'TMFast_1'; 'TMPref_1';'TMSlow_1'}; %add all trials here

%trialspeed=[1.4, 1.4, 1,4]; % double check this CURRENT VERSION DOESNT
%TAKE IN WALKING SPEED ANYWHERE

%Not Calculating step kinematics anymore in this version add this to run2
%Set up output matricies for SW(Step Width) and SL(Step Length): nconditions x 8
% output_sw_abs=zeros(12,8);output_sl_abs=zeros(12,8);
% step_width=cell(1,4); 
% step_length=cell(1,4);

for i=1:length(files)

    file=[char(files(i))];fl=char(files(i));
    
    %Load Data
    trc=load_trc(file,path); %reads marker kinematic data
    f=trc.freq(1); %Frequency of trc data 100hz for our system
    tdd=trc.time;%time scale for acceleration

    Analog=importdata([path, char(files(i)),'.anc'],'\t',11); %loads EMG data and GRF voltage data
    Analog=Analog.data;
    Names=importdata([path, char(files(i)),'.anc'],'\t',9); %Allows for picking out channels based on name
    
    trc.pos=trc.pos(:,:,:);

    Forces_temp=importdata([path, char(files(i)),'.forces'],'\t',5);
    Forces_temp=Forces_temp.data;
    Forces(i)={Forces_temp(:,:)};

    EMGchans=[width(Analog)-15:width(Analog)]; %EMG Channels are the last 16 for current setup. Could eventually edit to use system like below for finding ultrasonix column
    EMG(i)={Analog(:,EMGchans)}; 
    
     % Precision Stepping Trigger uses 
     in=('Ultrasonix');
     p=strcmpi(in,Names.textdata);
     [row,column] = find(p==1);
     
     if i==34 %Trial 34 is the only one that uses precision stepping trigger
            PrecisionTrigger = Analog(:,column);%Comment out for session 2
     end
    
    %Load Markers of interest (You have to edit this at the end too)
    S2_i=find(ismember(trc.mrk_names,'S2'));
    C7_i=find(ismember(trc.mrk_names,'C7'));
    rlr_i=find(ismember(trc.mrk_names,'R.Heel'));
    llr_i=find(ismember(trc.mrk_names,'L.Heel'));
    R_asis_i = find(ismember(trc.mrk_names,'R.ASIS'));
    L_asis_i = find(ismember(trc.mrk_names,'L.ASIS'));
    R_psis_i = find(ismember(trc.mrk_names,'R.PSIS'));
    L_psis_i = find(ismember(trc.mrk_names,'L.PSIS'));
    R_thigh_1_i = find(ismember(trc.mrk_names,'R.TH1'));
    R_thigh_2_i = find(ismember(trc.mrk_names,'R.TH2'));
    R_thigh_3_i = find(ismember(trc.mrk_names,'R.TH3'));
    R_knee_i = find(ismember(trc.mrk_names,'R.Knee'));
    R_shank_1_i = find(ismember(trc.mrk_names,'R.SH1'));
    R_shank_2_i = find(ismember(trc.mrk_names,'R.SH2'));
    R_shank_3_i = find(ismember(trc.mrk_names,'R.SH3'));
    R_shank_4_i = find(ismember(trc.mrk_names,'R.SH4'));
    R_ankle_i = find(ismember(trc.mrk_names,'R.Ankle'));
    R_mt_5_i = find(ismember(trc.mrk_names,'R.MT5'));
    R_mt_1_i = find(ismember(trc.mrk_names,'R.MT1'));
    L_thigh_1_i = find(ismember(trc.mrk_names,'L.TH1'));
    L_thigh_2_i = find(ismember(trc.mrk_names,'L.TH2'));
    L_thigh_3_i = find(ismember(trc.mrk_names,'L.TH3'));
    L_thigh_4_i = find(ismember(trc.mrk_names,'L.TH4'));
    L_knee_i = find(ismember(trc.mrk_names,'L.Knee'));
    L_shank_1_i = find(ismember(trc.mrk_names,'L.SH1'));
    L_shank_2_i = find(ismember(trc.mrk_names,'L.SH2'));
    L_shank_3_i = find(ismember(trc.mrk_names,'L.SH3'));
    L_ankle_i = find(ismember(trc.mrk_names,'L.Ankle'));
    L_mt_5_i = find(ismember(trc.mrk_names,'L.MT5'));
    L_mt_1_i = find(ismember(trc.mrk_names,'L.MT1'));
    R_acr_i = find(ismember(trc.mrk_names,'R.Acr'));
    L_acr_i = find(ismember(trc.mrk_names,'L.Acr'));
    T10_i = find(ismember(trc.mrk_names,'T10'));
    Clav_i = find(ismember(trc.mrk_names,'Clav'));
    Sternum_i = find(ismember(trc.mrk_names,'Sternum'));
    R_back_i = find(ismember(trc.mrk_names,'R.Back'));
    
    %Only trial 34 has Target 1 cleaned. Will through error for all other
    %trials if not seperated
     if i==34 %Precision Stepping 33 for YA01
         Target1_i = find(ismember(trc.mrk_names,'Target1'));
         Target1 = trc.pos(:,:,Target1_i);
         [b,a]=butter(4,16/f/2);
         targeta=filtfilt(b,a,Target1(:,:));
         StepTarget= {targeta};
     end
    
    
    %Pull in x,y,z data based on location
    C7=trc.pos(:,:,C7_i);
    rlr=trc.pos(:,:,rlr_i);
    llr=trc.pos(:,:,llr_i);
    sacrum=trc.pos(:,:,S2_i);
    R_asis = trc.pos(:,:,R_asis_i);
    L_asis = trc.pos(:,:,L_asis_i);
    R_psis = trc.pos(:,:,R_psis_i);
    L_psis = trc.pos(:,:,L_psis_i);
    R_thigh_1 = trc.pos(:,:,R_thigh_1_i);
    R_thigh_2 = trc.pos(:,:,R_thigh_2_i);
    R_thigh_3 = trc.pos(:,:,R_thigh_3_i);
    R_knee = trc.pos(:,:,R_knee_i);
    R_shank_1 = trc.pos(:,:,R_shank_1_i);
    R_shank_2 = trc.pos(:,:,R_shank_2_i);
    R_shank_3 = trc.pos(:,:,R_shank_3_i);
    R_shank_4 = trc.pos(:,:,R_shank_4_i);
    R_ankle = trc.pos(:,:,R_ankle_i);
    R_mt_5 = trc.pos(:,:,R_mt_5_i);
    R_mt_1 = trc.pos(:,:,R_mt_1_i);
    L_thigh_1 = trc.pos(:,:,L_thigh_1_i);
    L_thigh_2 = trc.pos(:,:,L_thigh_2_i);
    L_thigh_3 = trc.pos(:,:,L_thigh_3_i);
    L_thigh_4 = trc.pos(:,:,L_thigh_4_i);
    L_knee = trc.pos(:,:,L_knee_i);
    L_shank_1 = trc.pos(:,:,L_shank_1_i);
    L_shank_2 = trc.pos(:,:,L_shank_2_i);
    L_shank_3 = trc.pos(:,:,L_shank_3_i);
    L_ankle = trc.pos(:,:,L_ankle_i);
    L_mt_5 = trc.pos(:,:,L_mt_5_i);
    L_mt_1 = trc.pos(:,:,L_mt_1_i);
    R_acr = trc.pos(:,:,R_acr_i);
    L_acr = trc.pos(:,:,L_acr_i);
    T10 = trc.pos(:,:,T10_i);
    Clav = trc.pos(:,:,Clav_i);
    Sternum = trc.pos(:,:,Sternum_i);
    R_back = trc.pos(:,:,R_back_i);
    
    %Filter marker trajectories
    [b,a]=butter(4,16/f/2);
    
    rl=filtfilt(b,a,rlr(:,:));
    ll=filtfilt(b,a,llr(:,:));
    C7a=filtfilt(b,a,C7(:,:));
    s2a=filtfilt(b,a,sacrum(:,:));
    R_asisa = filtfilt(b,a,R_asis(:,:));
    L_asisa = filtfilt(b,a,L_asis(:,:));
    R_psisa = filtfilt(b,a,R_psis(:,:));
    L_psisa = filtfilt(b,a,L_psis(:,:));
    R_thigh_1a = filtfilt(b,a,R_thigh_1(:,:));
    R_thigh_2a = filtfilt(b,a,R_thigh_2(:,:));
    R_thigh_3a = filtfilt(b,a,R_thigh_3(:,:));
    R_kneea = filtfilt(b,a,R_knee(:,:));
    R_shank_1a = filtfilt(b,a,R_shank_1(:,:));
    R_shank_2a = filtfilt(b,a,R_shank_2(:,:));
    R_shank_3a = filtfilt(b,a,R_shank_3(:,:));
    R_shank_4a = filtfilt(b,a,R_shank_4(:,:));
    R_anklea = filtfilt(b,a,R_ankle(:,:));
    R_mt_5a = filtfilt(b,a,R_mt_5(:,:));
    R_mt_1a = filtfilt(b,a,R_mt_1(:,:));
    L_thigh_1a = filtfilt(b,a,L_thigh_1(:,:));
    L_thigh_2a = filtfilt(b,a,L_thigh_2(:,:));
    L_thigh_3a = filtfilt(b,a,L_thigh_3(:,:));
    L_thigh_4a = filtfilt(b,a,L_thigh_4(:,:));
    L_kneea = filtfilt(b,a,L_knee(:,:));
    L_shank_1a = filtfilt(b,a,L_shank_1(:,:));
    L_shank_2a = filtfilt(b,a,L_shank_2(:,:));
    L_shank_3a = filtfilt(b,a,L_shank_3(:,:));
    L_anklea = filtfilt(b,a,L_ankle(:,:));
    L_mt_5a = filtfilt(b,a,L_mt_5(:,:));
    L_mt_1a = filtfilt(b,a,L_mt_1(:,:));
    R_acra = filtfilt(b,a,R_acr(:,:));
    L_acra = filtfilt(b,a,L_acr(:,:));
    T10a = filtfilt(b,a,T10(:,:));
    Clava = filtfilt(b,a,Clav(:,:));
    Sternuma = filtfilt(b,a,Sternum(:,:));
    R_backa = filtfilt(b,a,R_back(:,:));
    
      
    %% Locate Heel Strikes
    %Depending on the direction participant is walking changes X or Y
    %direction of Heel/sacrum for heelstrikes
    if (i>5 && i<16) % Gait_init goes across width of lab
        [junk,rnn]=findpeaks(rl(:,2)-s2a(:,2));
        [junk,lnn]=findpeaks(ll(:,2)-s2a(:,2));   
        
    elseif (i>16 && i<33) %Overground walking flips direction back and forth so checks for that
        if s2a(1,1) < s2a(end,1)
            [junk,rnn]=findpeaks(rl(:,1)-s2a(:,1));
            [junk,lnn]=findpeaks(ll(:,1)-s2a(:,1)); 
        else
            [junk,rnn]=findpeaks(-rl(:,1)+s2a(:,1));
            [junk,lnn]=findpeaks(-ll(:,1)+s2a(:,1));
        end
        
    elseif i>33 %Normal treadmill walking
        [junk,rnn]=findpeaks(rl(:,1)-s2a(:,1));
        [junk,lnn]=findpeaks(ll(:,1)-s2a(:,1)); 
    end

    %% Locate Heel strikes
     if i>5 %Anything that isn't a fig 8 trial
        rn=[rnn zeros(length(rnn),1)];
        hsr=rnn(:,1);%heel-strike right

     
        ln=[lnn zeros(length(lnn),1)];
        hsl=lnn(:,1);%heel-strike left

     else %Finding heel strikes for fig 8s can't be done using normal method
         
         %heel strike position based off of z value once under minimum
         %value call good
         [junk,rnn]=findpeaks(-rl(:,3),'MinPeakProminence',0.01);
         var1 = -rl(:,3);
         ofst = rnn(end)+20;
         [standIn1,standIn2] = max(var1(ofst:end));
         rnn(end+1) = rnn(end)+20+standIn2;
         hsr=rnn(:,1);%heel-strike right
         
         
         [junk,lnn]=findpeaks(-ll(:,3),'MinPeakProminence',0.01);
         var1 = -ll(:,3);
         ofst = lnn(end)+20;
         [standIn1,standIn2] = max(var1(ofst:end));
         lnn(end+1) = lnn(end)+20+standIn2;
         hsl=lnn(:,1);%heel-strike right
         
     end
    
   
   %% Find Points
    lhs_temp=hsl(:,1);
    LHS(:,i)={lhs_temp};
    rhs_temp=hsr(:,1);
    RHS(:,i)={rhs_temp};

    Lheel(:,i)={ll};
    Rheel(:,i)={rl};
    Sacrum(:,i)={s2a};
    C7mark(:,i)={C7a};
    RAsis(:,i)= {R_asisa};
    LAsis(:,i)= {L_asisa};
    RPsis(:,i)= {R_psisa};
    LPsis(:,i)= {L_psisa};
    RTh1(:,i)= {R_thigh_1a};
    RTh2(:,i)= {R_thigh_2a};
    RTh3(:,i)= {R_thigh_3a};
    RKnee(:,i)= {R_kneea};
    RSh1(:,i)= {R_shank_1a};
    RSh2(:,i)= {R_shank_2a};
    RSh3(:,i)= {R_shank_3a};
    RSh4(:,i)= {R_shank_4a};
    RAnkle(:,i)= {R_anklea};
    RMT5(:,i)= {R_mt_5a};
    RMT1(:,i)= {R_mt_1a};
    LTh1(:,i)= {L_thigh_1a};
    LTh2(:,i)= {L_thigh_2a};
    LTh3(:,i)= {L_thigh_3a};
    LTh4(:,i)= {L_thigh_4a};
    LKnee(:,i)= {L_kneea};
    LSh1(:,i)= {L_shank_1a};
    LSh2(:,i)= {L_shank_2a};
    LSh3(:,i)= {L_shank_3a};
    LAnkle(:,i)= {L_anklea};
    LMT5(:,i)= {L_mt_5a};
    LMT1(:,i)= {L_mt_1a};
    RAcr(:,i)= {R_acra};
    LAcr(:,i)= {L_acra};
    T_10(:,i)= {T10a};
    Clavicle(:,i)= {Clava};
    Ster(:,i)= {Sternuma};
    RBack(:,i)= {R_back};

end

%% Compile All Output

for i=1:length(files)
     eval([['Output.' char(files(i))] '.LHS' '= LHS(:,i);'])
     eval([['Output.' char(files(i))] '.RHS' '= RHS(:,i);'])
    eval([['Output.' char(files(i))] '.Lheel' '= Lheel(:,i);'])
    eval([['Output.' char(files(i))] '.Rheel' '= Rheel(:,i);'])
    eval([['Output.' char(files(i))] '.Sacrum' '= Sacrum(:,i);'])
    eval([['Output.' char(files(i))] '.C7' '= C7mark(:,i);'])
    eval([['Output.' char(files(i))] '.RAsis' '= RAsis(:,i);'])
    eval([['Output.' char(files(i))] '.LAsis' '= LAsis(:,i);'])
    eval([['Output.' char(files(i))] '.RPsis' '= RPsis(:,i);'])
    eval([['Output.' char(files(i))] '.LPsis' '= LPsis(:,i);'])
    eval([['Output.' char(files(i))] '.RTh1' '= RTh1(:,i);'])
    eval([['Output.' char(files(i))] '.RTh2' '= RTh2(:,i);'])
    eval([['Output.' char(files(i))] '.RTh3' '= RTh3(:,i);'])
    eval([['Output.' char(files(i))] '.RKnee' '= RKnee(:,i);'])
    eval([['Output.' char(files(i))] '.RSh1' '= RSh1(:,i);'])
    eval([['Output.' char(files(i))] '.RSh2' '= RSh2(:,i);'])
    eval([['Output.' char(files(i))] '.RSh3' '= RSh3(:,i);'])
    eval([['Output.' char(files(i))] '.RSh4' '= RSh4(:,i);'])
    eval([['Output.' char(files(i))] '.RAnkle' '= RAnkle(:,i);'])
    eval([['Output.' char(files(i))] '.RMT5' '= RMT5(:,i);'])
    eval([['Output.' char(files(i))] '.RMT1' '= RMT1(:,i);'])
    eval([['Output.' char(files(i))] '.LTh1' '= LTh1(:,i);'])
    eval([['Output.' char(files(i))] '.LTh2' '= LTh2(:,i);'])
    eval([['Output.' char(files(i))] '.LTh3' '= LTh3(:,i);'])
    eval([['Output.' char(files(i))] '.LTh4' '= LTh4(:,i);'])
    eval([['Output.' char(files(i))] '.LKnee' '= LKnee(:,i);'])
    eval([['Output.' char(files(i))] '.LSh1' '= LSh1(:,i);'])
    eval([['Output.' char(files(i))] '.LSh2' '= LSh2(:,i);'])
    eval([['Output.' char(files(i))] '.LSh3' '= LSh3(:,i);'])
    eval([['Output.' char(files(i))] '.LAnkle' '= LAnkle(:,i);'])
    eval([['Output.' char(files(i))] '.LMT5' '= LMT5(:,i);'])
    eval([['Output.' char(files(i))] '.LMT5' '= LMT1(:,i);'])
    eval([['Output.' char(files(i))] '.RAcr' '= RAcr(:,i);'])
    eval([['Output.' char(files(i))] '.LAcr' '= LAcr(:,i);'])
    eval([['Output.' char(files(i))] '.T_10' '= T_10(:,i);'])
    eval([['Output.' char(files(i))] '.Clavicle' '= Clavicle(:,i);'])
    eval([['Output.' char(files(i))] '.Sternum' '= Ster(:,i);'])
    eval([['Output.' char(files(i))] '.RBack' '= RBack(:,i);'])
    
    
     if i==34 %Precision Stepping 33 for YA01
         eval([['Output.' char(files(i))] '.StepTarget' '= StepTarget;'])
         eval([['Output.' char(files(i))] '.PrecisionTrigger' '= PrecisionTrigger;'])
     end

    eval([['Output.' char(files(i))] '.EMG' '= EMG(i);'])
    eval([['Output.' char(files(i))] '.Forces' '= Forces(i);'])
end
%%
close all
clearvars -except Output
cd('D:\RESEARCH\Projects\Codes\Codes used in Per Motor Rep\SubjectOutput\') %Drag the output file into your Compiled folder and name it your subject number
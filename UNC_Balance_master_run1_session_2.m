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

%% Settings
sn='OANF10\Session2'; %Subject Number YA00\Session1
fs=100; %sampling frequency

% Define the path where your data is located:
path=char(['D:\RESEARCH\Projects\Data Collections\R21Repertoire\YA\clean\' sn '\Generated_files\']); %you will need to customized this.%E:\RESEARCH\Projects\Data Collections\R21Repertoire\YA\YA03\Session1

%Select Trials of Interest to Compile (list all as needed):

%  files = {'Fig8_1';'Fig8_2';'Fig8_3';'Fig8_4';'Fig8_5';'GaitInit_L_1';'GaitInit_L_2';'GaitInit_L_3';'GaitInit_L_4';'GaitInit_L_5';...
%     'GaitInit_R_1';'GaitInit_R_2';'GaitInit_R_3';'GaitInit_R_4';'GaitInit_R_5';'OGFast_1';'OGFast_2';'OGFast_3';'OGFast_4';'OGFast_5';'OGFast_6'...
%    ;'OGPref_1';'OGPref_2';'OGPref_3';'OGPref_4';'OGPref_5';'OGPref_6';'OGSlow_1';'OGSlow_2';'OGSlow_3';'OGSlow_4';'OGSlow_5';'OGSlow_6';...
%     'PrecisionStepping_1'; 'TMFast_1'; 'TMPref_1';'TMSlow_1'}; %add all trials here

%Session2
files = {'LFP_Left_1';'LFP_Right_1';'Norm_1'; 'Slip_1'; 'VR20_1'; 'VR35_1'}; %add all trials here

%Set EMG analog channels (open .anc and note the columns]:
%For Session 1
%EMGchans=[37:52];

%For Session 2
%EMGchans=[];


%Set walking speed for each/all trials:

%trialspeed=[1.4, 1.4, 1,4]; % double check this

%Set up output matricies for SW(Step Width) and SL(Step Length): nconditions x 8
% output_sw_abs=zeros(12,8);output_sl_abs=zeros(12,8);
% step_width=cell(1,4); 
% step_length=cell(1,4);

for i=1:length(files)
    %speed=trialspeed(i);
    file=[char(files(i))];fl=char(files(i));
    
    %Load Data
    trc=load_trc(file,path); %reads marker kinematic data
    f=trc.freq(1);
    tdd=trc.time;%time scale for acceleration

    Analog=importdata([path, char(files(i)),'.anc'],'\t',11); %loads EMG data and GRF voltage data
    Analog=Analog.data;
    Names=importdata([path, char(files(i)),'.anc'],'\t',9); %Allows for picking out channels based on name
    
    
    %Extract subregion of data
%     Analog=Analog(:,:);
    trc.pos=trc.pos(:,:,:);

    Forces_temp=importdata([path, char(files(i)),'.forces'],'\t',5);
    Forces_temp=Forces_temp.data;
    %Extract second minute for 2 min trials
%         tf=strncmp(files(i),'bf',2);
%         if tf==1
%             try
%               Analog=Analog(60001:120000,:);
%               Forces_temp=Forces_temp(60001:120000,:);
%               trc.pos=trc.pos(6001:12000,:,:);
%             end
%         end
        
    %Change file to 60s if collection ran long
%         tf=strcmp(files(i),'speed_11_1');
%         if tf==1
%           Analog=Analog(1:60000,:);
%           Forces_temp=Forces_temp(1:60000,:);
%           trc.pos=trc.pos(1:6000,:,:);
%         end

    Forces(i)={Forces_temp(:,:)};

    EMGchans=[width(Analog)-15:width(Analog)];
    EMG(i)={Analog(:,EMGchans)}; 
    
     % Precision Stepping Trigger
%      in=('Ultrasonix');
%      p=strcmpi(in,Names.textdata);
%      [row,column] = find(p==1);
%      
%      if i==34
%             PrecisionTrigger = Analog(:,column);%Comment out for session 2
%      end
    
    %Gait Initiation force plate data COMMENT OUT FOR SESSION 2
%      if i>5 && i<16
%          %step_side = Analog(:,column);      %Include this when we get to
%          %the trials that we start collecting when we say left/right
%          in=('F3X');
%          p=strcmpi(in,Names.textdata);
%          [row,column] = find(p==1);
%          GRFchans = [column:column+13];
%          GRF(i) = {Analog(:,GRFchans)};
%      end
    
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
    
    
%      if i==34 %Precision Stepping 33 for YA01
%          Target1_i = find(ismember(trc.mrk_names,'Target1'));
%          Target1 = trc.pos(:,:,Target1_i);
%          [b,a]=butter(4,16/f/2);
%          targeta=filtfilt(b,a,Target1(:,:));
%          StepTarget= {targeta};
%      end
    
    
    
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
    [junk,rnn]=findpeaks(rl(:,1)-s2a(:,1));
    [junk,lnn]=findpeaks(ll(:,1)-s2a(:,1)); 

    %% Locate 2nd point used for SW calculations
     rn=[rnn zeros(length(rnn),1)];
     hsr=rnn(:,1);%heel-strike right
%      for ii=1:length(rn)
%          if hsr(ii,1)+.75*f<length(rl)&&ii<length(rn)
%              rn(ii,2)=hsr(ii,1)+round(0.25.*(hsr(ii+1,1)-hsr(ii,1)));%find 25% of gc - ahead of initiation of terminal stance and heel rise
%              rn(ii,1)=hsr(ii,1)+round(0.12.*(hsr(ii+1,1)-hsr(ii,1)));%find 12% of gc
%          else
%              rn(ii,2)=rn(ii,1)+f*.25;
%              rn(ii,1)=rn(ii,1)+f*.12;
%              
%          end
%      end
     
     ln=[lnn zeros(length(lnn),1)];
     hsl=lnn(:,1);%heel-strike left
%      for ii=1:length(ln)
%          if ln(ii,1)+.75*f<length(ll)&&ii<length(ln)
%              ln(ii,2)=hsl(ii,1)+round(0.25.*(hsl(ii+1,1)-hsl(ii,1))); 
%              ln(ii,1)=hsl(ii,1)+round(0.12.*(hsl(ii+1,1)-hsl(ii,1))); 
%              
%          else
%              ln(ii,2)=hsl(ii,1)+f*.25;
%              ln(ii,2)=hsl(ii,1)+f*.12;
%          end
%      end
    
    %% Clean up 
%      drs=[];
%      
%      if rn(end,2)<120000 %should retain through 20 min (60*20*100)
%          drs=[rn];
%      else
%          drs=[rn(1:end-1,:)];
%      end
%      ddrs=drs(:,2)-drs(:,1);
%      [m,junk]=find(ddrs>f*.75);
%      drs(m,:)=[];
%      
%      dls=[];
%      
%      if ln(end,2)<120000 %should retain through 20 min (60*20*100)
%          dls=ln;
%      else
%          dls=[ln(1:end-1,:)];
%      end
%      ddls=dls(:,2)-dls(:,1);
%      [m,junk]=find(ddls>f*.75);
%      dls(m,:)=[];
    
   
%     figure
%     subplot(2,1,1)
%     plot((drs(:,1)+drs(:,2))./2,rv,'^r')
%     title(['Treadmill Velocity ' fl])
%     axis([0 18000 -2 0])
%     subplot(2,1,2)
%     plot((dls(:,1)+dls(:,2))./2,lv,'^k')
%     axis([0 18000 -2 0])

    
    %% Compute mean M-L position of heel markers for each stride over selected % of gait cycle
%     r=[];l=[];
%     drs=drs(1:end-1,:);dls=dls(1:end-1,:);
%     for j=1:length(drs)
%         x=mean(rl(drs(j,1):drs(j,2),1));
%         y=mean(rl(drs(j,1):drs(j,2),2));
%         z=mean(rl(drs(j,1):drs(j,2),3));
%         r.x(j,:)=x;
%         r.y(j,:)=y;
%         r.z(j,:)=z;
%     end
%     
%     for j=1:length(dls)
%         x=mean(ll(dls(j,1):dls(j,2),1));
%         y=mean(ll(dls(j,1):dls(j,2),2));
%         z=mean(ll(dls(j,1):dls(j,2),3));
%         l.x(j,:)=x;
%         l.y(j,:)=y;
%         l.z(j,:)=z;
%     end
       
    %% Reformat and compute step width
%     rmy=r.y;lmy=l.y;
%     
%     %Right and Left M-L variability
% 
%     drmy=diff(rmy);
%     mrwa=mean(drmy);srwa=std(drmy);
%     dlmy=diff(lmy);
%     mlwa=mean(dlmy);slwa=std(dlmy);
%        
%     ss=[drs (zeros(length(drs),1)) (1:length(drs))' rmy;dls ones(length(dls),1) (1:length(dls))' lmy];
%     ss=sortrows(ss,1);
%     
%     % Compute Step Width
%     sw=zeros(length(ss),1);
%     ds=diff(ss(:,3));
%     for ii=1:length(ss)-1
%         if ds(ii)>0
%             sw(ii,:)=ss(ii+1,5)-ss(ii,5);
%         elseif ds(ii)<0
%             sw(ii,:)=-(ss(ii+1,5)-ss(ii,5));
%         end
%     end
%    [em,junk]=find(sw==0);
%     if ~isempty(em)
%         sw(em)=[];
%     end
% 
%     step_width(1,i)={sw};
%     
%     mwa=mean(sw);%mean step width
%     dwa=std(sw);%standard deviation of step width over time
%  
%     output_sw_labels={'Mean Step Width';'Step Width StD';'Mean R.Heel';'R.Heel StD';'Mean L.Heel';'L.Heel StD';'# Right Steps';'# Left Steps'};
%     output_sw_abs(i,:)=[mwa dwa mrwa srwa mlwa slwa length(drs) length(dls)];
  
    
    %% Step length data
    
%     dtrl=rl(:,1)-(s2a(:,1));
%     dtll=ll(:,1)-s2a(:,1);
%     
%     rp=zeros(length(hsr)-1,1);lp=zeros(length(hsl)-1,1);
%     %Find 20% of GC for Right and Left strides
%     for ii=1:length(hsr)-1
%         rp(ii,1)=hsr(ii,1)+round(0.20.*(hsr(ii+1,1)-hsr(ii,1)));
%     end
%     for ii=1:length(hsl)-1
%         lp(ii,1)=hsl(ii,1)+round(0.20.*(hsl(ii+1,1)-hsl(ii,1)));
%     end
%     drl=zeros(length(rp)-1,1);dll=zeros(length(lp)-1,1);
%     %Compute Right Stride Length
%     for ii=1:length(rp)-1
%         drl(ii,1)=rl(rp(ii+1,1),1)-rl(rp(ii,1),1)+(speed.*(tdd(rp(ii+1,1))-tdd(rp(ii,1))));
%     end
%     %Compute Left Stride Length
%     for ii=1:length(lp)-1
%         dll(ii,1)=ll(lp(ii+1),1)-ll(lp(ii,1),1)+(speed.*(tdd(lp(ii+1),1)-tdd(lp(ii,1))));
%     end
%     
%     %Combine and sort for step length
%     pp=[rp ones(length(rp),1);lp zeros(length(lp),1)];pp=sortrows(pp,1);
%     dsl=zeros(length(pp)-1,1);
%     for ii=1:length(pp)-1
%         if pp(ii,2)==0
%             dsl(ii,1)=rl(pp(ii+1,1),1)-ll(pp(ii,1),1)+(speed.*(tdd(pp(ii+1,1))-tdd(pp(ii,1))));
%         elseif pp(ii,2)==1
%             dsl(ii,1)=ll(pp(ii+1,1),1)-rl(pp(ii,1),1)+(speed.*(tdd(pp(ii+1,1))-tdd(pp(ii,1))));
%         end
%         
%     end
% 
%     msa=mean(dsl);ssa=std(dsl); %Step mean and Variability
%     mra=mean(drl);sra=std(drl); %Rigth stride mean and variability
%     mla=mean(dll);sla=std(dll); %Left stride mean and variability
% 
%     step_length(1,i)={dsl};
% 
%     output_sl_labels={'Step Length Mean';'Step Length StD';'Right Stride Mean';'Right Stride StD';'Left Stride Mean';'Left Stride StD';'# Right Strides';'# Left Strides'};
%     output_sl_abs(i,:)=[msa ssa mra sra mla sla length(rp) length(lp)];   
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
    
    
%      if i==34 %Precision Stepping 33 for YA01
%          eval([['Output.' char(files(i))] '.StepTarget' '= StepTarget;'])
%          eval([['Output.' char(files(i))] '.PrecisionTrigger' '= PrecisionTrigger;'])
%      end
    % if i>5 && i<16
    %     eval([['Output.' char(files(i))] '.GRF' '= GRF(i);'])
     %end
    
    eval([['Output.' char(files(i))] '.EMG' '= EMG(i);'])
    eval([['Output.' char(files(i))] '.Forces' '= Forces(i);'])
   % eval([['Output.' char(files(i))] '.Stepwidth' '= step_width{1,i};'])
    %eval([['Output.' char(files(i))] '.Steplength' '= step_length{1,i};'])
end
%%
close all
clearvars -except Output
cd('D:\RESEARCH\Projects\Codes\Codes used in Per Motor Rep\SubjectOutput\') %Drag the output file into your Compiled folder and name it your subject number
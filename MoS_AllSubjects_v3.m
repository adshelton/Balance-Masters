%% Stand alone code for a participants Margin of Stability
%Current structure set for: Perturbation trials of Per Motor Rep 2nd
%session
%By: Andrew Shelton (Based on code by Brian Selgrade)
%EDITS:
%6/23/22 - ADS, Added in better comments

clear
clc
close all

%% Settings
fs=100; %sampling frequency
%All subjects interested in using
subjects = {'YA02', 'YA03', 'YA04', 'YA06', 'YA07', 'YA09', 'YA10', 'YA11', 'YA12', 'YA13','YA14','YA15','YA16','YA17','YA19','YA20','YA21','YA22','YA23','YA26'};
%file locations
p1 = 'D:\RESEARCH\Projects\Data Collections\R21Repertoire\YA\clean\';
p2 = '\Session2\Generated_files\';
%Select Trials of Interest to Compile:
files = {'Norm_1', 'VR20_1','VR35_1','Slip_1','LFP_Left_1','LFP_Right_1'};
% Initializing all of the directional MoS for each trial
Norm_APMosL = [];
Norm_MLMosL = [];
Norm_APMosR = [];
Norm_MLMosR = [];
VR20_APMosL = [];
VR20_MLMosL = [];
VR20_APMosR = [];
VR20_MLMosR = [];
VR35_APMosL = [];
VR35_MLMosL = [];
VR35_APMosR = [];
VR35_MLMosR = [];
SlipL_APMosR = [];
SlipL_MLMosR = [];
SlipR_APMosL = [];
SlipR_MLMosL = [];
LFPL_APMosL = [];
LFPL_MLMosL = [];
LFPR_APMosR = [];
LFPR_MLMosR = [];
for s = 1:length(subjects)
    path = strcat(p1, char(subjects(s)) ,p2); 
    clearvars -except fs path s subjects p1 p2 files MinMoSEarly_avg MinMoSEarly_std MinMoSLate_avg MinMoSLate_std MinMoSPost_avg MinMoSPost_std MinMoSBaseline_avg...
        MinMoSBaseline_std Norm_APMosL Norm_MLMosL Norm_APMosR Norm_MLMosR VR20_APMosL VR20_MLMosL VR20_APMosR VR20_MLMosR VR35_APMosL VR35_MLMosL VR35_APMosR...
        VR35_MLMosR SlipL_APMosR SlipL_MLMosR SlipR_APMosL SlipR_MLMosL LFPL_APMosL LFPL_MLMosL LFPR_APMosR LFPR_MLMosR Norm_APMos_avg Norm_MLMos_avg VR20_APMos_avg...
        VR20_MLMos_avg VR35_APMos_avg VR35_MLMos_avg Slip_APMos_first_avg Slip_MLMos_first_avg LFP_APMos_avg LFP_MLMos_avg Slip_APMos_p1_avg Slip_MLMos_p1_avg Slip_MLMos_p2_avg Slip_MLMos_p3_avg Slip_APMos_p2_avg Slip_APMos_p3_avg...
        Slip_APMos_p4_avg Slip_MLMos_p4_avg Slip_APMos_p5_avg Slip_MLMos_p5_avg Slip_APMos_p6_avg Slip_MLMos_p6_avg Slip_APMos_p7_avg Slip_MLMos_p7_avg Slip_APMos_p8_avg Slip_MLMos_p8_avg
    current = string(subjects(s));
    %Load in the heel strikes of interest for where perturbation occurs.
    %Done using slip_occur.m, LFP_occur_left.m, and LFP_occur_right.m
    load('D:\RESEARCH\Projects\Codes\Codes used in Per Motor Rep\HS\' + current + 'SlipHS');
    load('D:\RESEARCH\Projects\Codes\Codes used in Per Motor Rep\HS\' + current + 'LFP_Left_HS');
    load('D:\RESEARCH\Projects\Codes\Codes used in Per Motor Rep\HS\' + current + 'LFP_Right_HS');
    for i=1:length(files) %condition loop
        %next ~15 lines to load trc 
        clear LToeOff LP* LM* LS* LHS Lh* Ln* Ll* LA* RToeOff RP* RM* RS* RHS Rh* Rn* Rl* RA*
        file=char(files(i));

        %Load Data
        trc=load_trc(file,path); 
        f=trc.freq(1);
        tdd=trc.time;%time scale for acceleration
        trc.pos=trc.pos(:,:,:);

        %Load Markers of interest finds columns from trc file of each marker(You have to edit this at the end too)
        S2_i=find(ismember(trc.mrk_names,'S2'));
        RMT5_i=find(ismember(trc.mrk_names,'R.MT5'));
        LMT5_i=find(ismember(trc.mrk_names,'L.MT5'));
        RMT1_i=find(ismember(trc.mrk_names,'R.MT1'));
        LMT1_i=find(ismember(trc.mrk_names,'L.MT1'));
        rlr_i=find(ismember(trc.mrk_names,'R.Heel'));
        llr_i=find(ismember(trc.mrk_names,'L.Heel'));

        rlr=trc.pos(:,:,rlr_i);
        llr=trc.pos(:,:,llr_i);
        RMT1=trc.pos(:,:,RMT1_i);
        LMT1=trc.pos(:,:,LMT1_i);
        RMT5=trc.pos(:,:,RMT5_i);
        LMT5=trc.pos(:,:,LMT5_i);  
        sacrum=trc.pos(:,:,S2_i);

        %Filter marker trajectories
        [b,a]=butter(4,16/f/2);
        Rheel=filtfilt(b,a,rlr(:,:));
        Lheel=filtfilt(b,a,llr(:,:));
        s2a=filtfilt(b,a,sacrum(:,:));
        RMT1=filtfilt(b,a,RMT1(:,:));
        RMT5=filtfilt(b,a,RMT5(:,:));
        LMT1=filtfilt(b,a,LMT1(:,:));
        LMT5=filtfilt(b,a,LMT5(:,:));

    %% Calculating toe off and heel-strike
    g = 9.81; %gravity m/s^2
    VCoM = diff(s2a); % velocity of CoM

    RSacToe = s2a(:,1) - RMT1(:,1); %Distance from toe to sacrum. Used for toeoff
    LSacToe = s2a(:,1) - LMT1(:,1);

    RSacHeel = Rheel(:,1) - s2a(:,1); %Distance from sacrum to heel. Used for finding heelstrikes
    LSacHeel = Lheel(:,1) - s2a(:,1);

    [~, RToeOff] = findpeaks(RSacToe); %Calculate toe offs
    [~, LToeOff] = findpeaks(LSacToe);

    [~, RHS] = findpeaks(RSacHeel); %Calculate heel strikes. Outputs RHS as the frames where a heelstrike occurs
    [~, LHS] = findpeaks(LSacHeel);

    %% Calculate MoS
    % using method described by Young and Dingwell (2012)
    
    %XCoM
    %Following 2 for loops get leg length (distance from sacrum to heel in 3D space)
    for st = 1:length(RHS) % number of strides
        for ax = 1:3 % 3 axes
            diffR(st,ax) = Rheel(RHS(st),ax) - s2a(RHS(st),ax);
        end
        RL(st) = sqrt(diffR(st,1)^2 + diffR(st, 2)^2 + diffR(st, 3)^2);
    end

    for st = 1:length(LHS)
        for ax = 1:3
            diffL(st,ax) = Lheel(LHS(st),ax) - s2a(LHS(st),ax);
        end
        LL(st) = sqrt(diffL(st, 1)^2 + diffL(st, 2)^2 + diffL(st, 3)^2); 
    end

    meanLL = mean(LL); % left leg length
    meanRL = mean(RL); % right leg length

    Lavg = mean([meanRL, meanLL]); %leg length used (constant for each subject) 
    omega0 = sqrt(Lavg/g);

    XCoM = s2a(2:end,:) + VCoM*omega0; % extrapolated center of mass
    % BoS
    RMLBoS = RMT5(:,2);
    LMLBoS = LMT5(:,2);
    RAPBoS = RMT1(:,1);
    LAPBoS = LMT1(:,1);
    % MoS
    for k = 1:length(XCoM)
        LateralMoSRight(k,i) = 100*(XCoM(k,2) - RMLBoS(k+1)); % +1 due to velocity having 1 fewer point
        LateralMoSLeft(k,i) = 100*(LMLBoS(k+1) - XCoM(k,2));%100 converts MoS to cm
        AnteriorMoSRight(k,i) = 100*(RAPBoS(k+1) - XCoM(k,1));
        AnteriorMoSLeft(k,i) = 100*(LAPBoS(k+1) - XCoM(k,1));
    end

    
    %% Set value for each heelstrike
    %Section takes full MoS for each leg in ML and AP and takes out HS of
    %interest. All HS for continuous. Perturbation HS for discrete
    
    %Norm Trial
    % i corresponds with the trial currently on in the larger loop for the
    % subject
    j=0;
    if i==1
        for j = 1:length(LHS)
            Norm_APMosL(j) = AnteriorMoSLeft(LHS(j),i);
            Norm_MLMosL(j) = LateralMoSLeft(LHS(j),i);
        end
        j=0;
        for j = 1:length(RHS)
            Norm_APMosR(j) = AnteriorMoSRight(RHS(j),i);
            Norm_MLMosR(j) = LateralMoSRight(RHS(j),i);
        end
        j=0;
        %VR20
    elseif i==2
        for j = 1:length(LHS)
            VR20_APMosL(j) = AnteriorMoSLeft(LHS(j),i);
            VR20_MLMosL(j) = LateralMoSLeft(LHS(j),i);
        end
        j=0;
        for j = 1:length(RHS)
            VR20_APMosR(j) = AnteriorMoSRight(RHS(j),i);
            VR20_MLMosR(j) = LateralMoSRight(RHS(j),i);
        end
        j=0;
        %VR35
    elseif i==3
        for j = 1:length(LHS)
            VR35_APMosL(j) = AnteriorMoSLeft(LHS(j),i);
            VR35_MLMosL(j) = LateralMoSLeft(LHS(j),i);
        end
        j=0;
        for j = 1:length(RHS)
            VR35_APMosR(j) = AnteriorMoSRight(RHS(j),i);
            VR35_MLMosR(j) = LateralMoSRight(RHS(j),i);
        end
        j=0;
    % Slips - Left Slip want MoS for contraleral leg directly after   
    elseif i==4
        x=1; % x is for which number of the discrete perturbation this is
        for j = 1:length(LHS)
                standby = [];
            if any(SlipHS.Left == LHS(j))
                standby = find(RHS>LHS(j)); %Finding all contralateral that occur after the perturbation
                SlipL_APMosR_first(x) = AnteriorMoSRight(RHS(standby(1)),i); %Contraleteral recovery HS is the first that occurs in the standby variable
                SlipL_MLMosR_first(x) = LateralMoSRight(RHS(standby(1)),i);
                
                %Following lines are looking at the steps following the
                %discrete perturbation
                %P1 First left step on left side
                SlipL_APMosR_p1(x) = AnteriorMoSLeft(LHS(j+1),i);
                SlipL_MLMosR_p1(x) = LateralMoSLeft(LHS(j+1),i);
                %P2 Second Right step 
                SlipL_APMosR_p2(x) = AnteriorMoSRight(RHS(standby(2)),i);
                SlipL_MLMosR_p2(x) = LateralMoSRight(RHS(standby(2)),i);
                %P3 Second left step
                SlipL_APMosR_p3(x) = AnteriorMoSLeft(LHS(j+2),i);
                SlipL_MLMosR_p3(x) = LateralMoSLeft(LHS(j+2),i);
                %P4 Second Right step 
                SlipL_APMosR_p4(x) = AnteriorMoSRight(RHS(standby(3)),i);
                SlipL_MLMosR_p4(x) = LateralMoSRight(RHS(standby(3)),i);
                %P5 Second left step
                SlipL_APMosR_p5(x) = AnteriorMoSLeft(LHS(j+3),i);
                SlipL_MLMosR_p5(x) = LateralMoSLeft(LHS(j+3),i);
                %P6 Second Right step 
                SlipL_APMosR_p6(x) = AnteriorMoSRight(RHS(standby(4)),i);
                SlipL_MLMosR_p6(x) = LateralMoSRight(RHS(standby(4)),i);
                %P7 Second left step
                SlipL_APMosR_p7(x) = AnteriorMoSLeft(LHS(j+4),i);
                SlipL_MLMosR_p7(x) = LateralMoSLeft(LHS(j+4),i);
                %P8 Second Right step 
                SlipL_APMosR_p8(x) = AnteriorMoSRight(RHS(standby(5)),i);
                SlipL_MLMosR_p8(x) = LateralMoSRight(RHS(standby(5)),i);
                x=x+1;
                
            end
        end
        x=1;
        j=0;
        for j = 1:length(RHS)
                standby = [];
               %Same as above but for slips to right leg instead of left
            if any(SlipHS.Right == RHS(j))
                standby = find(LHS>RHS(j));
                SlipR_APMosL_first(x) = AnteriorMoSLeft(LHS(standby(1)),i);
                SlipR_MLMosL_first(x) = LateralMoSLeft(LHS(standby(1)),i);
                %P1 First left step on left side
                SlipR_APMosL_p1(x) = AnteriorMoSRight(RHS(j+1),i);
                SlipR_MLMosL_p1(x) = LateralMoSRight(RHS(j+1),i);
                %P2 Second Right step 
                SlipR_APMosL_p2(x) = AnteriorMoSLeft(LHS(standby(2)),i);
                SlipR_MLMosL_p2(x) = LateralMoSLeft(LHS(standby(2)),i);
                %P3 Second left step
                SlipR_APMosL_p3(x) = AnteriorMoSRight(RHS(j+2),i);
                SlipR_MLMosL_p3(x) = LateralMoSRight(RHS(j+2),i);
                %P4 Second Right step 
                SlipR_APMosL_p4(x) = AnteriorMoSLeft(LHS(standby(3)),i);
                SlipR_MLMosL_p4(x) = LateralMoSLeft(LHS(standby(3)),i);
                %P5 Second left step
                SlipR_APMosL_p5(x) = AnteriorMoSRight(RHS(j+3),i);
                SlipR_MLMosL_p5(x) = LateralMoSRight(RHS(j+3),i);
                %P6 Second Right step 
                SlipR_APMosL_p6(x) = AnteriorMoSLeft(LHS(standby(4)),i);
                SlipR_MLMosL_p6(x) = LateralMoSLeft(LHS(standby(4)),i);
                %P7 Second left step
                SlipR_APMosL_p7(x) = AnteriorMoSRight(RHS(j+4),i);
                SlipR_MLMosL_p7(x) = LateralMoSRight(RHS(j+4),i);
                %P8 Second Right step 
                SlipR_APMosL_p8(x) = AnteriorMoSLeft(LHS(standby(5)),i);
                SlipR_MLMosL_p8(x) = LateralMoSLeft(LHS(standby(5)),i);
                
                x=x+1;
            end
        end 
        x=1;
        j=0;
    % LFP_Left    
    elseif i==5
        for j = 1:length(LHS)
            %LFP_occur codes already give HS and recovery step is on same
            %side so need for standby variable
         if any(LFP_Left_HS.Swing == LHS(j))
               LFPL_APMosL(x) = AnteriorMoSLeft(LHS(j+1),i); % LFP_Left gives HS at the start of GC involving left pull
               LFPL_MLMosL(x) = LateralMoSLeft(LHS(j+1),i); %Want MoS at heel strike directly after LFP
               x=x+1;
         end
        end
        x=1;
         j=0;
    %LFP_Right     
    elseif i==6
        x=1;
        for j = 1:length(RHS)
         if any(LFP_Right_HS.Swing == RHS(j))
               LFPR_APMosR(x) = AnteriorMoSRight(RHS(j+1),i); % LFP_Left gives HS at the start of GC involving left pull
               LFPR_MLMosR(x) = LateralMoSRight(RHS(j+1),i); %Want MoS at heel strike directly after LFP
               x=x+1;
         end
        end
    end
    
   

    end
    %% Outcomes
    %Chance that there are zero MoS values (Really only for if there are not 5 slips or waist pulls. Or if we don't get the full lsit of recovery steps for slips/waist pulls at the end before data collection stops)
    Norm_APMosL(Norm_APMosL == 0) = NaN;
    Norm_APMosR(Norm_APMosR == 0) = NaN;
    Norm_MLMosL(Norm_MLMosL == 0) = NaN;
    Norm_MLMosR(Norm_MLMosR == 0) = NaN;
    
    %Averages together MoS for left and right side to give one directional
    %value per particpant per trial
    Norm_APMos_avg(s) = mean([Norm_APMosL, Norm_APMosR]);
    Norm_MLMos_avg(s) = mean([Norm_MLMosL, Norm_MLMosR]);
    
    VR20_APMosL(VR20_APMosL == 0) = NaN;
    VR20_APMosR(VR20_APMosR == 0) = NaN;
    VR20_MLMosL(VR20_MLMosL == 0) = NaN;
    VR20_MLMosR(VR20_MLMosR == 0) = NaN;
    
    VR20_APMos_avg(s) = mean([VR20_APMosL, VR20_APMosR]);
    VR20_MLMos_avg(s) = mean([VR20_MLMosL, VR20_MLMosR]);
    
    VR35_APMosL(VR35_APMosL == 0) = NaN;
    VR35_APMosR(VR35_APMosR == 0) = NaN;
    VR35_MLMosL(VR35_MLMosL == 0) = NaN;
    VR35_MLMosR(VR35_MLMosR == 0) = NaN;
    
    VR35_APMos_avg(s) = mean([VR35_APMosL, VR35_APMosR]);
    VR35_MLMos_avg(s) = mean([VR35_MLMosL, VR35_MLMosR]);
    
    SlipL_APMosR_first(SlipL_APMosR_first == 0) = NaN;
    SlipR_APMosL_first(SlipR_APMosL_first == 0) = NaN;
    SlipL_MLMosR_first(SlipL_MLMosR_first == 0) = NaN;
    SlipR_MLMosL_first(SlipR_MLMosL_first == 0) = NaN;
    
    Slip_APMos_first_avg(s) = mean([SlipL_APMosR_first, SlipR_APMosL_first]);
    Slip_MLMos_first_avg(s) = mean([SlipL_MLMosR_first, SlipR_MLMosL_first]);
    
    SlipL_APMosR_p1(SlipL_APMosR_p1 == 0) = NaN;
    SlipR_APMosL_p1(SlipR_APMosL_p1 == 0) = NaN;
    SlipL_MLMosR_p1(SlipL_MLMosR_p1 == 0) = NaN;
    SlipR_MLMosL_p1(SlipR_MLMosL_p1 == 0) = NaN;
    
    Slip_APMos_p1_avg(s) = mean([SlipL_APMosR_p1, SlipR_APMosL_p1]);
    Slip_MLMos_p1_avg(s) = mean([SlipL_MLMosR_p1, SlipR_MLMosL_p1]);
    
    SlipL_APMosR_p2(SlipL_APMosR_p2 == 0) = NaN;
    SlipR_APMosL_p2(SlipR_APMosL_p2 == 0) = NaN;
    SlipL_MLMosR_p2(SlipL_MLMosR_p2 == 0) = NaN;
    SlipR_MLMosL_p2(SlipR_MLMosL_p2 == 0) = NaN;
    
    Slip_APMos_p2_avg(s) = mean([SlipL_APMosR_p2, SlipR_APMosL_p2]);
    Slip_MLMos_p2_avg(s) = mean([SlipL_MLMosR_p2, SlipR_MLMosL_p2]);
    
    SlipL_APMosR_p3(SlipL_APMosR_p3 == 0) = NaN;
    SlipR_APMosL_p3(SlipR_APMosL_p3 == 0) = NaN;
    SlipL_MLMosR_p3(SlipL_MLMosR_p3 == 0) = NaN;
    SlipR_MLMosL_p3(SlipR_MLMosL_p3 == 0) = NaN;
    
    Slip_APMos_p3_avg(s) = mean([SlipL_APMosR_p3, SlipR_APMosL_p3]);
    Slip_MLMos_p3_avg(s) = mean([SlipL_MLMosR_p3, SlipR_MLMosL_p3]);
    
    SlipL_APMosR_p4(SlipL_APMosR_p4 == 0) = NaN;
    SlipR_APMosL_p4(SlipR_APMosL_p4 == 0) = NaN;
    SlipL_MLMosR_p4(SlipL_MLMosR_p4 == 0) = NaN;
    SlipR_MLMosL_p4(SlipR_MLMosL_p4 == 0) = NaN;
    
    Slip_APMos_p4_avg(s) = mean([SlipL_APMosR_p4, SlipR_APMosL_p4]);
    Slip_MLMos_p4_avg(s) = mean([SlipL_MLMosR_p4, SlipR_MLMosL_p4]);
   
    SlipL_APMosR_p5(SlipL_APMosR_p5 == 0) = NaN;
    SlipR_APMosL_p5(SlipR_APMosL_p5 == 0) = NaN;
    SlipL_MLMosR_p5(SlipL_MLMosR_p5 == 0) = NaN;
    SlipR_MLMosL_p5(SlipR_MLMosL_p5 == 0) = NaN;
    
    Slip_APMos_p5_avg(s) = mean([SlipL_APMosR_p5, SlipR_APMosL_p5]);
    Slip_MLMos_p5_avg(s) = mean([SlipL_MLMosR_p5, SlipR_MLMosL_p5]);
    
    SlipL_APMosR_p6(SlipL_APMosR_p6 == 0) = NaN;
    SlipR_APMosL_p6(SlipR_APMosL_p6 == 0) = NaN;
    SlipL_MLMosR_p6(SlipL_MLMosR_p6 == 0) = NaN;
    SlipR_MLMosL_p6(SlipR_MLMosL_p6 == 0) = NaN;
    
    Slip_APMos_p6_avg(s) = mean([SlipL_APMosR_p6, SlipR_APMosL_p6]);
    Slip_MLMos_p6_avg(s) = mean([SlipL_MLMosR_p6, SlipR_MLMosL_p6]);
    
    SlipL_APMosR_p7(SlipL_APMosR_p7 == 0) = NaN;
    SlipR_APMosL_p7(SlipR_APMosL_p7 == 0) = NaN;
    SlipL_MLMosR_p7(SlipL_MLMosR_p7 == 0) = NaN;
    SlipR_MLMosL_p7(SlipR_MLMosL_p7 == 0) = NaN;
    
    Slip_APMos_p7_avg(s) = mean([SlipL_APMosR_p7, SlipR_APMosL_p7]);
    Slip_MLMos_p7_avg(s) = mean([SlipL_MLMosR_p7, SlipR_MLMosL_p7]);
    
    SlipL_APMosR_p8(SlipL_APMosR_p8 == 0) = NaN;
    SlipR_APMosL_p8(SlipR_APMosL_p8 == 0) = NaN;
    SlipL_MLMosR_p8(SlipL_MLMosR_p8 == 0) = NaN;
    SlipR_MLMosL_p8(SlipR_MLMosL_p8 == 0) = NaN;
    
    Slip_APMos_p8_avg(s) = mean([SlipL_APMosR_p8, SlipR_APMosL_p8]);
    Slip_MLMos_p8_avg(s) = mean([SlipL_MLMosR_p8, SlipR_MLMosL_p8]);
    
    LFPL_APMosL(LFPL_APMosL == 0) = NaN;
    LFPR_APMosR(LFPR_APMosR == 0) = NaN;
    LFPL_MLMosL(LFPL_MLMosL == 0) = NaN;
    LFPR_MLMosR(LFPR_MLMosR == 0) = NaN;
    
    LFP_APMos_avg(s) = mean([LFPL_APMosL, LFPR_APMosR]);
    LFP_MLMos_avg(s) = mean([LFPL_MLMosL, LFPR_MLMosR]);
    
    %Clear variables to prep for next run through of code
    clearvars Norm_APMosL Norm_MLMosL Norm_APMosR Norm_MLMosR VR20_APMosL VR20_MLMosL VR20_APMosR VR20_MLMosR VR35_APMosL VR35_MLMosL VR35_APMosR...
        VR35_MLMosR SlipL_APMosR_first SlipL_MLMosR_first SlipR_APMosL_first SlipR_MLMosL_first...
        SlipL_APMosR_p1 SlipL_MLMosR_p1 SlipR_APMosL_p1 SlipR_MLMosL_p1...
        SlipL_APMosR_p2 SlipL_MLMosR_p2 SlipR_APMosL_p2 SlipR_MLMosL_p2...
        SlipL_APMosR_p3 SlipL_MLMosR_p3 SlipR_APMosL_p3 SlipR_MLMosL_p3...
        SlipL_APMosR_p4 SlipL_MLMosR_p4 SlipR_APMosL_p4 SlipR_MLMosL_p4...
        SlipL_APMosR_p5 SlipL_MLMosR_p5 SlipR_APMosL_p5 SlipR_MLMosL_p5...
        SlipL_APMosR_p6 SlipL_MLMosR_p6 SlipR_APMosL_p6 SlipR_MLMosL_p6...
        SlipL_APMosR_p7 SlipL_MLMosR_p7 SlipR_APMosL_p7 SlipR_MLMosL_p7...
        SlipL_APMosR_p8 SlipL_MLMosR_p8 SlipR_APMosL_p8 SlipR_MLMosL_p8...
        LFPL_APMosL LFPL_MLMosL LFPR_APMosR LFPR_MLMosR
end

% Clearing out some more variables. Loops already run through so not
% entirely necessary
clearvars -except Norm_APMos_avg Norm_MLMos_avg VR20_APMos_avg VR20_MLMos_avg VR35_APMos_avg VR35_MLMos_avg Slip_APMos_first_avg Slip_MLMos_first_avg...
    LFP_APMos_avg LFP_MLMos_avg Slip_APMos_p1_avg Slip_MLMos_p1_avg Slip_MLMos_p2_avg Slip_MLMos_p3_avg Slip_APMos_p2_avg Slip_APMos_p3_avg...
    Slip_APMos_p4_avg Slip_MLMos_p4_avg Slip_APMos_p5_avg Slip_MLMos_p5_avg Slip_APMos_p6_avg Slip_MLMos_p6_avg Slip_APMos_p7_avg Slip_MLMos_p7_avg Slip_APMos_p8_avg Slip_MLMos_p8_avg

%% Plotting Boxplots
%removed VR20
n = 20;
x = [repmat(1,n,1) repmat(2,n,1) repmat(3,n,1) repmat(4,n,1) ];
 
figure(1)
boxplot([Norm_APMos_avg', VR35_APMos_avg', Slip_APMos_first_avg', LFP_APMos_avg'],'Labels',{'Norm', 'VR35','Slip','Waist Pull'})
hold on;
title('AP MoS')
ylabel('Margin ofStability (cm)')
data1 = [Norm_APMos_avg' VR35_APMos_avg' Slip_APMos_first_avg' LFP_APMos_avg'];
scatter(x(:),data1(:),'filled','MarkerFaceAlpha', 0.6','jitter','on','jitterAmount',0.025);
hold on;

figure(2)
boxplot([Norm_MLMos_avg', VR35_MLMos_avg', Slip_MLMos_first_avg', LFP_MLMos_avg'],'Labels',{'Norm', 'VR35','Slip','Waist Pull'})
hold on;
title('ML MoS')
ylabel('Margin ofStability (cm)')
data2 = [Norm_MLMos_avg' VR35_MLMos_avg' Slip_MLMos_first_avg' LFP_MLMos_avg'];
scatter(x(:),data2(:),'filled','MarkerFaceAlpha', 0.6','jitter','on','jitterAmount',0.025);
hold on;

%% Stats
% ANOVA
% Just do it in R

%% WBAM Slip Progression code
% Used for prelimanry WBAM data to look at how WBAM Range changes as
% participants get more used to the slip perturbation
%Andrew Shelton
%3/5/2022

%% Setup

clear;
clc;
close all;

subj = {'YA02','YA03','YA04','YA06', 'YA07', 'YA09', 'YA10', 'YA11', 'YA12','YA13', 'YA14', 'YA15', 'YA16', 'YA17','YA19','YA20', 'YA21','YA22','YA23','YA26'};

%% Code for output files
% Take in all WBAM/break down into sections/do area between the curves
for i= 1:length(subj)
    %Load necessary files
    %Subject .mat files
    current = string(subj(i));
    load('D:\RESEARCH\Projects\Codes\Codes used in Per Motor Rep\SubjectOutput\' + current + 'Session2');
    load('D:\RESEARCH\Projects\Codes\Codes used in Per Motor Rep\HS\' + current + 'SlipHS');
    load('D:\RESEARCH\Projects\Codes\Codes used in Per Motor Rep\HS\' + current + 'LFP_Left_HS');
    load('D:\RESEARCH\Projects\Codes\Codes used in Per Motor Rep\HS\' + current + 'LFP_Right_HS');
    
    %All WBAM Files
    wbam_norm = load('D:\RESEARCH\Projects\Codes\WBAM\PerMotorRepBOSSubData\' + current + '\Norm_1_WBAM');
    wbam_slip = load('D:\RESEARCH\Projects\Codes\WBAM\PerMotorRepBOSSubData\' + current + '\Slip_1_WBAM');
    wbam_vr20 = load('D:\RESEARCH\Projects\Codes\WBAM\PerMotorRepBOSSubData\' + current + '\VR20_1_WBAM');
    wbam_vr35 = load('D:\RESEARCH\Projects\Codes\WBAM\PerMotorRepBOSSubData\' + current + '\VR35_1_WBAM');
    wbam_left = load('D:\RESEARCH\Projects\Codes\WBAM\PerMotorRepBOSSubData\' + current + '\LFP_Left_1_WBAM');
    wbam_right = load('D:\RESEARCH\Projects\Codes\WBAM\PerMotorRepBOSSubData\' + current + '\LFP_Right_1_WBAM');
    
    
    %---SLIP--- (Left foot slips)
    lhs = Output.Slip_1.LHS{1,1};
    stride_avg_sag = [];
    stride_avg_trans = [];
    stride_avg_front = [];
    start = 0;
    last = 0;
    values_stride_sag = [];
    values_stride_trans = [];
    values_stride_front = [];

    percent_GC_sag = [];
    percent_GC_trans = [];
    percent_GC_front = [];
    sub_avg_sag = 0;
    sub_avg_trans = 0;
    sub_avg_front = 0;
    GC_sag = [];
    GC_trans = [];
    GC_front = [];
    counter = 0;
    for j = 1:length(lhs)-1
        start = lhs(j);
        last = lhs(j+1);
        if any(SlipHS.Left == start)
            counter = counter + 1;
            %Pull out WBAM values
            values_stride_sag = wbam_slip.ans(start:last,91);
            values_stride_trans = wbam_slip.ans(start:last,90);
            values_stride_front = wbam_slip.ans(start:last,89);
            values_stride_front = values_stride_front.NWBAM_R .*-1;
            %Range for each stride
            stride_avg_sag_l(i,counter) = max(values_stride_sag.NWBAM_F)-min(values_stride_sag.NWBAM_F);
            stride_avg_trans_l(i,counter) = max(values_stride_trans.NWBAM_U)-min(values_stride_trans.NWBAM_U);
            stride_avg_front_l(i, counter) = max(values_stride_front)-min(values_stride_front);
            
            
        end 
        
    end
   

    
    %---SLIP--- (Right foot splips)
    rhs = Output.Slip_1.RHS{1,1};
    stride_avg_sag = [];
    stride_avg_trans = [];
    stride_avg_front = [];
    start = 0;
    last = 0;
    values_stride_sag = [];
    values_stride_trans = [];
    values_stride_front = [];
    
    percent_GC_sag = [];
    percent_GC_trans = [];
    percent_GC_front = [];
    sub_avg_sag = 0;
    sub_avg_trans = 0;
    sub_avg_front = 0;
    GC_sag = [];
    GC_trans = [];
    GC_front = [];
    counter = 0;
    for j = 1:length(rhs)-1
        start = rhs(j);
        last = rhs(j+1);
        if any(SlipHS.Right == start)
            counter = counter+1
            %Pull out WBAM values
            values_stride_sag = wbam_slip.ans(start:last,91);
            values_stride_trans = wbam_slip.ans(start:last,90);
            values_stride_front = wbam_slip.ans(start:last,89);
            values_stride_front = values_stride_front.NWBAM_R .*-1;
            %Range for each stride
            stride_avg_sag_r(i,counter) = max(values_stride_sag.NWBAM_F)-min(values_stride_sag.NWBAM_F);
            stride_avg_trans_r(i,counter) = max(values_stride_trans.NWBAM_U)-min(values_stride_trans.NWBAM_U);
            stride_avg_front_r(i,counter) = max(values_stride_front)-min(values_stride_front);
            
            
        end   
    end
    
    
end

%% PLOTS
stride_avg_sag_l(stride_avg_sag_l==0) = nan;
stride_avg_trans_l(stride_avg_trans_l==0) = nan;
stride_avg_front_l(stride_avg_front_l==0) = nan;

stride_avg_sag_r(stride_avg_sag_r==0) = nan;
stride_avg_trans_r(stride_avg_trans_r==0) = nan;
stride_avg_front_r(stride_avg_front_r==0) = nan;
   
figure(1)
boxplot(stride_avg_sag_l)
title('Left Sagittal')
ylim([0.05 0.20])
figure(2)
boxplot(stride_avg_trans_l)
title('Left Transverse')
ylim([0.010 0.040])
figure(3)
boxplot(stride_avg_front_l)
title("Left Frontal")
ylim([0.01 0.08])
figure(4)
boxplot(stride_avg_sag_r)
title("Right Sagittal")
ylim([0.05 0.20])
figure(5)
boxplot(stride_avg_trans_r)
title("Right Transverse")
ylim([0.010 0.040])
figure(6)
boxplot(stride_avg_front_r)            
title("Right Frontal")
ylim([0.01 0.08])








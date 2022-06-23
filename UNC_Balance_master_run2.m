%% Master analayis routine 2 for walking data and balance analyses
% Uses kinematics to define gait cycle events
% CURRENTLY WIP by ADS as of 6/23/22 to update to current Run1 outputs
%Edits needed: Calculate SW, calculate SL(maybe?), what is a hodograph?,
%cleanup divergence exponent section to just be one run instead of needing
%also run3
clear
clc

% Input: 
%       Individial subject data, compiled in Run1
% Optional Output: 
%       Step kinematics and gait variability (sk)
%       Short- and long-term divergence exponents (le)
%       Margins of stability (mos)
%       Hodographs (hg)


fs=100;

% Setup output (0=no, 1=yes)
sk=1;
le=1;
mos=1;
hg=1;

%Trials for your subjects
conds={'Fixed_NoPert_1_clean','Fixed_Pert_1_clean','SP_NoPert_1_clean','SP_Pert_1_clean'};
% conds={'VR_0_125','VR_20_125','VR_35_125','VR_50_125'};
filenames = uigetfile('*.*','Select All Files','multiselect','on'); %Opens up file explorer and select all SUbject outputs that you want to run

%% Runs through Outputs of choice 

for s=1:length(filenames)
   
    load(filenames{1,s});
%     load(filenames)
    for c=1:length(conds)
        cond=[char(conds(c))];
        cond=['Output.' cond];

        data=eval(cond);

%         if sk==1; %Step Kinematics    
%             SW=data.Stepwidth*100;
%             SL=data.Steplength*100;
%             SWavg(s,c)=nanmean(SW);
%             SLavg(s,c)=nanmean(SL);
%             SWsd(s,c)=nanstd(SW);
%             SLsd(s,c)=nanstd(SL);
%             S2_sdx(s,c)=std(data.Sacrum(:,1));
%             S2_sdy(s,c)=std(data.Sacrum(:,2));
%             S2_sdz(s,c)=std(data.Sacrum(:,3));
% %             S2V_sdx(s,c)=std(fs*[diff(data.Sacrum(:,1));0]);
% %             S2V_sdy(s,c)=std(fs*[diff(data.Sacrum(:,2));0]);
% %             S2V_sdz(s,c)=std(fs*[diff(data.Sacrum(:,3));0]);
%         end

%         if mos==1;
%             RHS{s,c}=data.RHS{1,1}(:,1);
%             LHS{s,c}=data.LHS{1,1}(:,1);
% %             data.Forces{1,1}(data.Forces{1,1}==0)=NaN;
% %             CoPx1{s,c}=data.Forces{1,1}(:,6);
% %             CoPy1{s,c}=data.Forces{1,1}(:,5);
% %             CoPx2{s,c}=data.Forces{1,1}(:,13);
% %             CoPy2{s,c}=data.Forces{1,1}(:,12);
%             S2{s,c}=data.Sacrum(:,:);
%             Rheel{s,c}=data.Rheel(:,:);
%             Lheel{s,c}=data.Lheel(:,:);
%         end
        
        
%         if hg==1;
% %             RHS=data.RHS{1,1}(:,1);
% %             LHS=data.LHS{1,1}(:,1);
%             for p=1:3;
%                 S2=data.Sacrum(:,:);
%                 C7=data.C7(:,:);
% %                 S2v(:,p)=fs*[diff(S2(:,p)); 0];
%             end
%             
%             nstrides=length(RHS);
%             for i=1:nstrides-1
% %                 S2v_x_temp=S2v(RHS(i):RHS(i+1),1);
% %                 S2v_y_temp=S2v(RHS(i):RHS(i+1),2);
% %                 S2v_z_temp=S2v(RHS(i):RHS(i+1),3);
%                 
% %                 t=[1/length(S2v_x_temp):1/length(S2v_x_temp):1];
%                 tt=[1/1000:1/1000:1];
% %                 S2v_x_strides(:,i)=interp1(t,S2v_x_temp,tt);
% %                 S2v_y_strides(:,i)=interp1(t,S2v_y_temp,tt);
% %                 S2v_z_strides(:,i)=interp1(t,S2v_z_temp,tt);
%             end
%             
% %             S2v_x(:,s,c)=nanmean(S2v_x_strides,2);
% %             S2v_y(:,s,c)=nanmean(S2v_y_strides,2);
% %             S2v_z(:,s,c)=nanmean(S2v_z_strides,2);
%         end
        
        if le==1
%             ts=[0:(length(Lheel{s,c}(:,1))-1)]/fs;
            
            %for p=1:3
%                 Lheelv{s,c}{:,p}=fs*[diff(Lheel{s,c}{:,p})];
%                 Rheelv{s,c}{:,p}=fs*[diff(Rheel{s,c}{:,p})];
%                 S2v{s,c}{:,p}=fs*[diff(S2{s,c}{:,p})];
%                 C7v{s,c}{:,p}=fs*[diff(C7{s,c}{:,p})];
            %end
            x = data.RHS{1,1}(:,1);
            if s==1 && c==1
                logicalIndexes = x < 30000 & x > 12000;
            else
                logicalIndexes = x < 33000 & x > 15000;
            end
            RHS{s,c}=x(logicalIndexes);
            %S2{s,c}=data.Sacrum(:,:);
            m=4; 
            meanperiod=round(mean(diff(RHS{s,c})));
            maxiter=1100;
            tao=round(meanperiod/4); 
            
%             state=[S2{s,c}(:,:),S2v{s,c}(:,:)];
            state=cell2mat(data.Sacrum(:,:));
            if s==1 && c==1
                state = state(12000:30000,:);
            else
                state = state(15000:33000,:);
            end
%             state=[C7v{s,c}(:,:)]; %3D trunk
% %              state=[C7v{s,c}(:,1:3)-S2v{s,c}(:,1:3)]; % State parameter for dynamic stability
%             state=[S2v{s,c}(:,3)];
            
            fprintf('\n');
            [filenames(s),conds(c)]
            d_all{s,c} = lyarosenstein(state,m,tao,meanperiod,maxiter);
            mn_pd(s,c)=meanperiod;
            tau_all(s,c)=tao;

        end
        
    end
    
    
    
end


%% Compute group-average values



% for c=1:length(conds)
%     if sk==1;
%         GroupOutput.SWavg(c)=mean(SWavg(:,c));GroupOutput.SWsd(c)=std(SWavg(:,c));%/sqrt(length(filenames));
%         GroupOutput.SLavg(c)=mean(SLavg(:,c));GroupOutput.SLsd(c)=std(SLavg(:,c));%/sqrt(length(filenames));
%         GroupOutput.SWVavg(c)=mean(SWsd(:,c));GroupOutput.SWVsd(c)=std(SWsd(:,c));%/sqrt(length(filenames));
%         GroupOutput.SLVavg(c)=mean(SLsd(:,c));GroupOutput.SLVsd(c)=std(SLsd(:,c));%/sqrt(length(filenames));
%         GroupOutput.S2_sdx_avg(c)=mean(S2_sdx(:,c));GroupOutput.S2_sdx_sd(c)=std(S2_sdx(:,c));%/sqrt(length(filenames));
%         GroupOutput.S2_sdy_avg(c)=mean(S2_sdy(:,c));GroupOutput.S2_sdy_sd(c)=std(S2_sdy(:,c));%/sqrt(length(filenames));
%         GroupOutput.S2_sdz_avg(c)=mean(S2_sdz(:,c));GroupOutput.S2_sdz_sd(c)=std(S2_sdz(:,c));%/sqrt(length(filenames));
%         GroupOutput.S2V_sdx_avg(c)=mean(S2V_sdx(:,c));GroupOutput.S2V_sdx_sd(c)=std(S2V_sdx(:,c));%/sqrt(length(filenames));
%         GroupOutput.S2V_sdy_avg(c)=mean(S2V_sdy(:,c));GroupOutput.S2V_sdy_sd(c)=std(S2V_sdy(:,c));%/sqrt(length(filenames));
%         GroupOutput.S2V_sdz_avg(c)=mean(S2V_sdz(:,c));GroupOutput.S2V_sdz_sd(c)=std(S2V_sdz(:,c));%/sqrt(length(filenames));
%     end
%     
%     if hg==1;
%         for i=1:1000;
%             S2v_x_avg(i,c)=nanmean(S2v_x(i,:,c));
%             S2v_y_avg(i,c)=nanmean(S2v_y(i,:,c));
%             S2v_z_avg(i,c)=nanmean(S2v_z(i,:,c));
%         end
%     end
% end

%% Plot group-average values
% Cond=[1:5]
% Cond2=[Cond]+.1;
% if sk==1
%     figure(1)
%     plot(SWsd(:,1),100*S2_sdy(:,1),'ko'), hold on
%     plot(SWsd(:,2:4),100*S2_sdy(:,2:4),'ro')
%     
%     figure(2)
%     plot(SWsd(:,2:4)-SWsd(:,1),100*S2_sdy(:,2:4)-S2_sdy(:,1),'ro')
% %    [RHO,PVAL] =corr([(SWsd(:,2)-SWsd(:,1));(SWsd(:,3)-SWsd(:,1));(SWsd(:,2)-SWsd(:,1))],[(S2_sdy(:,2)-S2_sdy(:,1));(S2_sdy(:,3)-S2_sdy(:,1));(S2_sdy(:,2)-S2_sdy(:,1))])
%     [RHO_swpost,PVAL_swpost] = corr((SWsd(:,4)-SWsd(:,1)),(100*S2_sdy(:,4)-S2_sdy(:,1)))
%    
%     figure(3)
%     subplot(121)
%     plot(SWsd(:,1),SWsd(:,4)-SWsd(:,1),'ko'), hold on
%     plot(SWsd(:,1),SWsd(:,3)-SWsd(:,1),'ro'), hold on
%     plot(SWsd(:,1),SWsd(:,2)-SWsd(:,1),'bo'), hold on
%     
%     [RHO_sw(1),PVAL_sw(1)]=corr(SWsd(:,1),(SWsd(:,2)-SWsd(:,1)));
%     [RHO_sw(2),PVAL_sw(2)]=corr(SWsd(:,1),(SWsd(:,3)-SWsd(:,1)));
%     [RHO_sw(3),PVAL_sw(3)]=corr(SWsd(:,1),(SWsd(:,4)-SWsd(:,1)))
%     
%     subplot(122)
%     plot(100*S2_sdy(:,1),100*S2_sdy(:,4)-100*S2_sdy(:,1),'ko')
%     
%     [RHO_post(1),PVAL_post(1)]=corr(S2_sdy(:,1),(S2_sdy(:,2)-S2_sdy(:,1)));
%     [RHO_post(2),PVAL_post(2)]=corr(S2_sdy(:,1),(S2_sdy(:,3)-S2_sdy(:,1)));
%     [RHO_post(3),PVAL_post(3)]=corr(S2_sdy(:,1),(S2_sdy(:,4)-S2_sdy(:,1)))
    
%      figure
%     subplot(221)
%     plot(Cond,GroupOutput.SLavg(1:5),'ko'), hold on
%     errorbar(Cond,GroupOutput.SLavg(1:5),GroupOutput.SLsd(1:5),'k')
%     ylim([0 80])
%     ylabel('Stride Length(cm)')
%     plot(Cond2,GroupOutput.SLavg(6:10),'ro'), hold on
%     errorbar(Cond2,GroupOutput.SLavg(6:10),GroupOutput.SLsd(6:10),'r')
%     ylim([0 80])
% 
%     subplot(222)
%     plot(Cond,GroupOutput.SWavg(1:5),'ko'), hold on
%     errorbar(Cond,GroupOutput.SWavg(1:5),GroupOutput.SWsd(1:5),'k')
%     ylim([0 20])
%     ylabel('Stride Width (cm)')
%     plot(Cond2,GroupOutput.SWavg(6:10),'ro'), hold on
%     errorbar(Cond2,GroupOutput.SWavg(6:10),GroupOutput.SWsd(6:10),'r')
%     ylim([0 20])
% 
%     subplot(223)
%     plot(Cond,GroupOutput.SLVavg(1:5),'ko'), hold on
%     errorbar(Cond,GroupOutput.SLVavg(1:5),GroupOutput.SLVsd(1:5),'k')
%     ylim([0 4])
%     ylabel('Stride Length Variability (cm)')
%     plot(Cond2,GroupOutput.SLVavg(6:10),'ro'), hold on
%     errorbar(Cond2,GroupOutput.SLVavg(6:10),GroupOutput.SLVsd(6:10),'r')
%     ylim([0 4])
% 
%     Cond=[1:5];
%     subplot(224)
%     plot(Cond,GroupOutput.SWVavg(1:5),'ko'), hold on
%     errorbar(Cond,GroupOutput.SWVavg(1:5),GroupOutput.SWVsd(1:5),'k')
%     ylim([0 4])
%     ylabel('Stride Width Variability (cm)')
%     plot(Cond2,GroupOutput.SWVavg(6:10),'ro'), hold on
%     errorbar(Cond2,GroupOutput.SWVavg(6:10),GroupOutput.SWVsd(6:10),'r')
%     ylim([0 4])
%     
%     figure
%     subplot(131)
%     plot(Cond,GroupOutput.S2_sdx_avg(1:5),'ko'), hold on
%     errorbar(Cond,GroupOutput.S2_sdx_avg(1:5),GroupOutput.S2_sdx_sd(1:5),'k')
%     ylabel('S2_sdx(cm)')
%     plot(Cond2,GroupOutput.S2_sdx_avg(6:10),'ro'), hold on
%     errorbar(Cond2,GroupOutput.S2_sdx_avg(6:10),GroupOutput.S2_sdx_sd(6:10),'r')
%     ylim([0 0.08])
% 
%     subplot(132)
%     plot(Cond,GroupOutput.S2_sdy_avg(1:5),'ko'), hold on
%     errorbar(Cond,GroupOutput.S2_sdy_avg(1:5),GroupOutput.S2_sdy_sd(1:5),'k')
%     ylabel('S2_sdy(cm)')
%     plot(Cond2,GroupOutput.S2_sdy_avg(6:10),'ro'), hold on
%     errorbar(Cond2,GroupOutput.S2_sdy_avg(6:10),GroupOutput.S2_sdy_sd(6:10),'r')
%     ylim([0 0.08])
% 
%     subplot(133)
%     plot(Cond,GroupOutput.S2_sdz_avg(1:5),'ko'), hold on
%     errorbar(Cond,GroupOutput.S2_sdz_avg(1:5),GroupOutput.S2_sdz_sd(1:5),'k')
%     ylabel('S2_sdz(cm)')
%     plot(Cond2,GroupOutput.S2_sdz_avg(6:10),'ro'), hold on
%     errorbar(Cond2,GroupOutput.S2_sdz_avg(6:10),GroupOutput.S2_sdz_sd(6:10),'r')
%     ylim([0 0.08])
% 
% 
% figure
%     subplot(131)
%     plot(Cond,GroupOutput.S2V_sdx_avg(1:5),'ko'), hold on
%     errorbar(Cond,GroupOutput.S2V_sdx_avg(1:5),GroupOutput.S2V_sdx_sd(1:5),'k')
%     ylabel('S2V_sdx(cm)')
%     plot(Cond2,GroupOutput.S2V_sdx_avg(6:10),'ro'), hold on
%     errorbar(Cond2,GroupOutput.S2V_sdx_avg(6:10),GroupOutput.S2V_sdx_sd(6:10),'r')
%     ylim([0 0.3])
% 
%     subplot(132)
%     plot(Cond,GroupOutput.S2V_sdy_avg(1:5),'ko'), hold on
%     errorbar(Cond,GroupOutput.S2V_sdy_avg(1:5),GroupOutput.S2V_sdy_sd(1:5),'k')
%     ylabel('S2V_sdy(cm)')
%     plot(Cond2,GroupOutput.S2V_sdy_avg(6:10),'ro'), hold on
%     errorbar(Cond2,GroupOutput.S2V_sdy_avg(6:10),GroupOutput.S2V_sdy_sd(6:10),'r')
%     ylim([0 0.3])
% 
%     subplot(133)
%     plot(Cond,GroupOutput.S2V_sdz_avg(1:5),'ko'), hold on
%     errorbar(Cond,GroupOutput.S2V_sdz_avg(1:5),GroupOutput.S2V_sdz_sd(1:5),'k')
%     ylabel('S2V_sdz(cm)')
%     plot(Cond2,GroupOutput.S2V_sdz_avg(6:10),'ro'), hold on
%     errorbar(Cond2,GroupOutput.S2V_sdz_avg(6:10),GroupOutput.S2V_sdz_sd(6:10),'r')
%     ylim([0 0.3])
    
%end

% if hg==1;
%     figure
%     subplot(251), plot(S2v_x_avg(15:end,1),S2v_z_avg(15:end,1)),axis([-.3 .3 -.3 .3])
%     subplot(252), plot(S2v_x_avg(15:end,2),S2v_z_avg(15:end,2)),axis([-.3 .3 -.3 .3])
%     subplot(253), plot(S2v_x_avg(15:end,3),S2v_z_avg(15:end,3)),axis([-.3 .3 -.3 .3])
%     subplot(254), plot(S2v_x_avg(15:end,4),S2v_z_avg(15:end,4)),axis([-.3 .3 -.3 .3])
%     subplot(255), plot(S2v_x_avg(15:end,5),S2v_z_avg(15:end,5)),axis([-.3 .3 -.3 .3])
%     
%     subplot(256), plot(S2v_x_avg(15:end,6),S2v_z_avg(15:end,6)),axis([-.3 .3 -.3 .3])
%     subplot(257), plot(S2v_x_avg(15:end,7),S2v_z_avg(15:end,7)),axis([-.3 .3 -.3 .3])
%     subplot(258), plot(S2v_x_avg(15:end,8),S2v_z_avg(15:end,8)),axis([-.3 .3 -.3 .3])
%     subplot(259), plot(S2v_x_avg(15:end,9),S2v_z_avg(15:end,9)),axis([-.3 .3 -.3 .3])
%     subplot(2,5,10), plot(S2v_x_avg(15:end,10),S2v_z_avg(15:end,10)),axis([-.3 .3 -.3 .3])
%     
%     figure
%     subplot(251), plot(S2v_y_avg(15:end,1),S2v_x_avg(15:end,1)),axis([-.3 .3 -.3 .3])
%     subplot(252), plot(S2v_y_avg(15:end,2),S2v_x_avg(15:end,2)),axis([-.3 .3 -.3 .3])
%     subplot(253), plot(S2v_y_avg(15:end,3),S2v_x_avg(15:end,3)),axis([-.3 .3 -.3 .3])
%     subplot(254), plot(S2v_y_avg(15:end,4),S2v_x_avg(15:end,4)),axis([-.3 .3 -.3 .3])
%     subplot(255), plot(S2v_y_avg(15:end,5),S2v_x_avg(15:end,5)),axis([-.3 .3 -.3 .3])
%     
%     subplot(256), plot(S2v_y_avg(15:end,6),S2v_x_avg(15:end,6)),axis([-.3 .3 -.3 .3])
%     subplot(257), plot(S2v_y_avg(15:end,7),S2v_x_avg(15:end,7)),axis([-.3 .3 -.3 .3])
%     subplot(258), plot(S2v_y_avg(15:end,8),S2v_x_avg(15:end,8)),axis([-.3 .3 -.3 .3])
%     subplot(259), plot(S2v_y_avg(15:end,9),S2v_x_avg(15:end,9)),axis([-.3 .3 -.3 .3])
%     subplot(2,5,10), plot(S2v_y_avg(15:end,10),S2v_x_avg(15:end,10)),axis([-.3 .3 -.3 .3])
%     
%     
%     
%     
% end

%%
% if mos==1
% %     figure
% %     plot([0.001:.001:30],CoPx1{2,10}(1:30000)), hold on
% %     plot([0.001:.001:30],CoPx2{2,10}(1:30000))
% %     plot([0.01:.01:30],S2{2,10}(1:3000)*1000)
% 
% 
%     figure
%     plot(Rheel{1,1}(:,1)*1000,Rheel{1,1}(:,2)*1000,'r'), hold on
%     plot(CoPy1{1,1}(:,1),CoPx1{1,1}(:,1),'k')
%     plot(Lheel{1,1}(:,1)*1000,Lheel{1,1}(:,2)*1000,'b')
%     plot(CoPy2{1,1}(:,1),CoPx2{1,1}(:,1),'g')
% 
% 
%     plot(S2{1,1}(:,1)*1000,S2{1,1}(:,2)*1000,'c'), hold on
% 
% 
% end

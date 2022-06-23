% Post processing of divergence calculations
% LLE Calculation
%clear all; close all;

% load LLE_S2_3d_pandv.mat
% load LLE_S2_3d_v.mat
% load LLE_S2V_z.mat
% load LLE_old_pilot.mat

% for i=1:3
% if i==1;load LLE_S2V_x.mat;end
% if i==2;load LLE_S2V_y.mat;end
% if i==3;load LLE_S2V_z.mat;end


try;tao/meanperiod,end
try; comment{:}, end
try; note{:}, end
sz=[1,1];%size(filenames)

cols=get(gca,'ColorOrder');
lysty={'-',':','--','-.','-',':','--','-',':','--','-.','-',':','--'};
for fnum=1:11%:length(filenames);
%     figure(fnum); clf(fnum);
    for cnum=1:4 %Number of conditions
    
    d=d_all{fnum,cnum};%(:,fnum,cnum);
    d=d';
%     lh(cnum)=line([1:length(d)],d,'LineStyle',lysty{cnum});%'Color',cols(cnum,:)
%     xlabel('Sample #'); ylabel('Divergence');
%     title(filenames(fnum));

    tlinear_s=[1:100];
    tlinear_l=[400:1000];
    F_s = polyfit(tlinear_s(:),d(tlinear_s),1);
%     line(tlinear_s,polyval(F_s,tlinear_s),'LineWidth',2,'LineStyle',lysty{cnum});%'Color',cols(cnum,:)
    lle_s(fnum,cnum) = F_s(1)*mn_pd(fnum); % normalize as Dingwell did
    
    
    polydata = polyval(F_s,tlinear_s(:));
    sstot = sum((d(tlinear_s) - mean(d(tlinear_s))).^2);
    ssres = sum((d(tlinear_s) - polydata).^2);
    rsquared_s(fnum,cnum) = 1 - (ssres / sstot);
    
    
    
    F_l = polyfit(tlinear_l(:),d(tlinear_l),1);
%     line(tlinear_l,polyval(F_l,tlinear_l),'LineWidth',2,'LineStyle',lysty{cnum});%'Color',cols(cnum,:)
    lle_l(fnum,cnum) = F_l(1)*mn_pd(fnum);

    polydata = polyval(F_l,tlinear_l(:));
    sstot = sum((d(tlinear_l) - mean(d(tlinear_l))).^2);
    ssres = sum((d(tlinear_l) - polydata).^2);
    rsquared_l(fnum,cnum) = 1 - (ssres / sstot);
    
    
    
    
    % 
%     %gait cycle specific Lyapunov Exponents
%         figure(31);
%         plot(d(1:5,:).'); hold on;
%         plot(d(6:10,:).',':'); hold off
%         xlabel('Sample #'); ylabel('log(Divergence)');

        % [xrange,yrange]=ginput(2);
        % tlinear=round(xrange(1)):round(xrange(2));
%         tlinear=1:100;%400:1000; % 20:100;%
%         for k=1:10;
%             F(k,:) = polyfit(tlinear,d(k,tlinear),1);
%             line(tlinear,polyval(F(k,:),tlinear),'Color','r','LineWidth',2);
%             lle(k,fnum,cnum) = F(k,1)*meanperiod;
%         end
       
    end
end

%%
% close all

for c=1:4 % number of conditions
    
Group.lle_s_avg(:,c)=mean(lle_s(:,c));
Group.lle_s_sd(:,c)=std(lle_s(:,c));

Group.lle_l_avg(:,c)=mean(lle_l(:,c));
Group.lle_l_sd(:,c)=std(lle_l(:,c));

end

figure
Cond=[1:4];
subplot(121)
plot(Cond,Group.lle_s_avg(1:4),'ko'), hold on
% errorbar(Cond,Group.lle_s_avg(1:5),Group.lle_s_sd(1:5),'k')
%axis([0 7 0 1.2])
% plot(Cond+.1,Group.lle_s_avg(6:10),'ro'), hold on
% errorbar(Cond+.1,Group.lle_s_avg(6:10),Group.lle_s_sd(6:10),'r')
plot(3,Group.lle_s_avg(4),'ro'), hold on

%ylim([0 1.2])
subplot(122)
plot(Cond,Group.lle_l_avg(1:4),'ko'), hold on
errorbar(Cond,Group.lle_l_avg(1:4),Group.lle_l_sd(1:4),'k')
%ylim([0 .15])
% plot(Cond+0.1,Group.lle_l_avg(6:10),'ro'), hold on
% errorbar(Cond+0.1,Group.lle_l_avg(6:10),Group.lle_l_sd(6:10),'r')
%ylim([0 .15])
% end
%% Compute and Plot Average divergence curves
% figure(3)
% for i=1:length(d_all(:,:,1))
%     d_mean(i,1)=(mean(d_all(i,:,1)));
%     d_se(i,1)=(std(d_all(i,:,1)))/sqrt(10);
%     
%     d_mean(i,2)=(mean(d_all(i,:,2)));
%     d_se(i,2)=(std(d_all(i,:,2)))/sqrt(10);
%     
%     d_mean(i,3)=(mean(d_all(i,:,3)));
%     d_se(i,3)=(std(d_all(i,:,3)))/sqrt(10);
%     
%     d_mean(i,4)=(mean(d_all(i,:,4)));
%     d_se(i,4)=(std(d_all(i,:,4)))/sqrt(10);
%     
%     d_mean(i,5)=(mean(d_all(i,:,5)));
%     d_se(i,5)=(std(d_all(i,:,5)))/sqrt(10);
%     
%     d_mean(i,6)=(mean(d_all(i,:,6)));
%     d_se(i,6)=(std(d_all(i,:,6)))/sqrt(10);
% 
% end
% 
% 
% 
% plot(d_mean)
% subplot(211)
% plot(d_young_mean), hold on
% legend('normal', 'visual')
% plot(d_young_mean-d_young_se,'--')
% plot(d_young_mean+d_young_se,'--')
% xlim([0 1000])
% ylabel('< ln[dj(t)] >')
% title('YOUNG')
% 
% % axis([0 1000 -5 0])
% subplot(212)
% plot(d_old_mean),hold on
% plot(d_old_mean-d_old_se,'--')
% plot(d_old_mean+d_old_se,'--')
% xlim([0 1000])
% % axis([0 1000 -5 0])
% ylabel('< ln[dj(t)] >')
% title('OLD')
% xlabel('time')
% 
% figure(4)
% subplot(211)
% plot(exp(d_young_mean));
% legend('normal', 'visual')
% % plot(d_old_mean-d_old_se,'--')
% % plot(d_old_mean+d_old_se,'--')
% title('YOUNG')
% xlim([0 1000])
% ylabel('< dj (t) >')
% 
% 
% subplot(212)
% plot(exp(d_old_mean));
% xlim([0 1000])
% xlabel('time')
% ylabel('< dj (t) >')
% title('OLD')
% % 
% % 
% % 

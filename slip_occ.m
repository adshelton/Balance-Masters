%% WHEN DOES SLIP OCCUR?
%Pretty slef explanatory title. Give HS that triggers the slip perturbation
%to occur
%Andrew Shelton 11/03/21
%
clear
clc

%Load File
load('D:\RESEARCH\Projects\Codes\Codes used in Per Motor Rep\SubjectOutput\OA01Session2') %Need to run UNC_Balance_Master_run1 for this first

%% Left Leg
%Call in position
x = Output.Slip_1.Lheel{1,1}(:,1);

%Velocity calculation
delta_t = 0.01;
t_l = length(x);
v=diff(x)./delta_t;
t = 0:0.01:(t_l*0.01);

%Acceleration Calcs
a=diff(v)./delta_t; % accelerations at times ta; a vector of the length less than t, x by 2

%Find when the slips are occuring
[var, loc] = findpeaks(-a); %Dececleration so flip the acceleration so find peaks function works
slips_pre = loc(find((var<10) & (var>3))); %Slips are  at "6 m/s2" but marker acceleration doesn't always show that 

%Determine if the slip is really for this leg (sometime heels both accelerate at the same speed depending on how the participant is recovering to slip)
slips = [];
for i = 1:length(slips_pre)
   if (Output.Slip_1.Lheel{1,1}(slips_pre(i),3) < Output.Slip_1.Rheel{1,1}(slips_pre(i),3)) %Which foot is on the ground at acceleration
       slips(end+1) = slips_pre(i);
   end
end

%Heel strike for step where slip occurs
strike = Output.Slip_1.LHS{1,1}; %All footstrikes. Need to run UNC_Balance_Master_run1 for this first
j =1;
distance = 100;
%Point of acceleration peak doesn't always line exactly up with HS.
%Slightly slower so need this code to pick up which HS frame exactly
%triggered it. Helps to line up with later MoS code
strike_slip = [];
for i = 1:length(strike)
    distance_1 = slips(j) - strike(i);
    if((distance > distance_1) && distance_1>0)
        distance = distance_1;
    end
    if(distance_1 < 0)
       strike_slip(end+1) = strike(i-1);
       distance = 100;
       j = j + 1;
    end
    if j > length(slips)
        break
    end
end

strike_slip_L = strike_slip

%% Right Leg
%Same as above but for slips that are occuring to Right leg
%Call in position
x = Output.Slip_1.Rheel{1,1}(:,1);

%Velocity calculation
delta_t = 0.01;
t_l = length(x);
v=diff(x)./delta_t;
t = 0:0.01:(t_l*0.01);

%Acceleration Calcs
a=diff(v)./delta_t; % accelerations at times ta; a vector of the length less than t, x by 2

%Find when the slips are occuring
[var, loc] = findpeaks(-a);
slips_pre = loc(find((var<9) & (var>3.5)));

%Determine if the slip is really for this leg
slips = [];
for i = 1:length(slips_pre)
   if (Output.Slip_1.Rheel{1,1}(slips_pre(i),3) < Output.Slip_1.Lheel{1,1}(slips_pre(i),3))
       slips(end+1) = slips_pre(i);
   end
end

%Heel strike for step where slip occurs
strike = Output.Slip_1.RHS{1,1};
j =1;
distance = 100;
strike_slip = [];
for i = 1:length(strike)
    distance_1 = slips(j) - strike(i);
    if((distance > distance_1) && distance_1>0)
        distance = distance_1;
    end
    if(distance_1 < 0)
       strike_slip(end+1) = strike(i-1);
       distance = 100;
       j = j + 1;
    end
    if j > length(slips)
        break
    end
end

strike_slip_R = strike_slip

%Make Matlab structure
SlipHS.Left = strike_slip_L;
SlipHS.Right = strike_slip_R;

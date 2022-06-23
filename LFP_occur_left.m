%% When does LFP occur?
% Used to determine the stride where the perturbation occurs for lateral
% waist pulls. Code for pulls to left leg at toe off. Flip variables for
% right leg
% By: Andrew Shelton
% 6/23/22: Added comments

%% Initialize
clear
clc

%Load File for current subject
load('D:\RESEARCH\Projects\Codes\Codes used in Per Motor Rep\SubjectOutput\OA10Session2')

j = 1;   %Number of heels strikes found
distance = 100;     %Distance variable for determing what heel strike goes along with the perturbation
strike_pull_swing = []; %Heel strikes that start stride where participants are pulled towards their swing leg at toe off (Left leg when pulled to right)
strike_pull_stance = []; %Heel strike that starts stride where the leg is in stance during perturbation (Right leg when pulled to left)

%% Left Leg
%Frame where participant has taken a wider recovery step (using Cortex marker X,Y,Z graph just some point from recovery step is in stance)
pulls = [544;1079;1719;2365;3008];

%Heel strike for step where slip occurs
strike = Output.LFP_Left_1.LHS{1,1};    %Heel strikes for side that participant is pulled towards
stance_strike = Output.LFP_Left_1.RHS{1,1}; %Opposite side foot strikes

%Loop to find all strikes relevant to perturbations (probably a cleaner way to do this, but made sense at the time)
for i = 1:length(strike) %Cycle through all heel strikes
    distance_1 = pulls(j) - strike(i);  %New distance of frames between when a pull occured and when a strike occured
    if((distance > distance_1) && distance_1>0) %(To be completely honest need to double check this if statement. Don't think it does anything)
        distance = distance_1;
    end
    if(distance_1 < 0) % When distance first goes negative it means that that is the heelstrike right after the recovery step
       strike_pull_swing(end+1) = strike(i-2); % want the heel strike that starts swing phase (right after recovery step so need to subtract two for start)
       strike_pull_stance(end+1) = stance_strike(find((stance_strike > strike(i-2)) & (stance_strike < strike(i-1)))); %Finding corresponding stance heelstrike
       distance = 100; %Reset distance parameter
       j = j + 1; %Go to next pull
    end
    if j > length(pulls)        %Break out of for loop when all perturbation heel strikes have been found
        break
    end
end

% Set variables to matlab structure
LFP_Left_HS.Swing = strike_pull_swing
LFP_Left_HS.Stance = strike_pull_stance
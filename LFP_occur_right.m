%% When does LFP occur
%Same as LFP_occur_left except flipped directions. Left code better
%commented.
%By: Andrew Shelton


clear
clc
%Load File
load('D:\RESEARCH\Projects\Codes\Codes used in Per Motor Rep\SubjectOutput\OA10Session2') % Participant Output 

%% Right Leg
%When pulls occur
pulls = [751;1726;2696;3772;4731]; % Manual entry of where the pulls occur (see image for where these points come from)


%Heel strike for step where slip occurs
strike = Output.LFP_Right_1.RHS{1,1}; %Right foot heel strikes (Switch if pulling to left)
stance_strike = Output.LFP_Right_1.LHS{1,1}; % Left foot heel strikes (Switch if pulling to left)
j =1;
distance = 100;
strike_pull_swing = [];
strike_pull_stance = [];
for i = 1:length(strike)
    distance_1 = pulls(j) - strike(i);  %First section of loop determines what heel strike occurs closest to that frame
    if((distance > distance_1) && distance_1>0)   
        distance = distance_1;
    end
    if(distance_1 < 0)
       strike_pull_swing(end+1) = strike(i-2); % want the heel strike that starts swing phase so the last right leg heel strike before perturbation on right side
       strike_pull_stance(end+1) = stance_strike(find((stance_strike > strike(i-2)) & (stance_strike < strike(i-1)))); %Gives heel strike for left leg right before perturbation on right side occurs
       distance = 100;
       j = j + 1;
    end
    if j > length(pulls) %Once all heel strikes are found exits out of code
        break
    end
end

%Builds output structure - Will need to then be saved in your desired
%location
LFP_Right_HS.Swing = strike_pull_swing 
LFP_Right_HS.Stance = strike_pull_stance
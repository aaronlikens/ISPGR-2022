
%% TITLE: MFDFATutorial.m 
% DATE: June, 2022
% AUTHOR: Anaelle E. Charles, MS
% EMAIL: bmchnonan@unomaha.edu

% UPDATED: June 28th, 2022

% DESCRIPTION:
% Section I: - Tutorial Example
%            - Contain step by step tutorial to run MFDFA Analysis.

% Section II: - Full analysis
%             - Load all 3 subjects
%             - Run the analysis for all subjects through the loops
%             - Create summary table

% Copyright 2022,  Anaelle E. Charles

% Redistribution and use of this script, with or without
% modification, is permitted provided this copyright notice,
% the original authors name and the following disclaimer remains.

%% Required Toolboxes:
% Statistics and Machine Learning Toolbox
% Signal Processing Toolbox
% Image Processing Toolbox
%% ------------------ DOWNLOAD MFDFA FOLDER -------------------------------
    %Option 1: Clone Nonan lib and Fork tutorial materials on personal desktop
    %Futher instruction [insert here]  https://github.com/Nonlinear-Analysis-Core
    %Option 2: Download- MFDFA_Workshop folder on your desktop
 
addpath('MFDFA')

my_directory = uigetdir(matlabroot, 'select StrideInterval folder');

%% Section 1: Step-by-step Tutorial Example
    
%%---------------------Step 1:Load Stride Interval-------------------------

% Source of data: Raffalt et al., 2021
%Load data: 
load S206_selfPaced_StrideIntervals.mat

time= SI(:,1);
StrideIntervals= SI(:,2);

%% ------------Step 2: Visual inspection before analysis:------------------
    %   Plot and inspect your data:
        %   It's important that you look at your graph before analysis.
        %   Make sure that you don't see any abnormal large peaks or other
        %   abnormalities.

    %   What to look for?
        %   At this point of the pre-processing, your data should look
        %   noise-like time series (refer to slide 25-26 of the
        %   power-point for best practices).

figure(1);
plot(StrideIntervals, 'r-')
title('Walking Stride intervals time series', 'FontSize', 16)
xlabel('Time (s)');
ylabel('Stride Interval (s)');
% Enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
disp('Press any key to continue!')

%% ------------Step 3: Set MFDFA parameters for Analysis-------------------
%MFDFA function parameters: Details can be found in the MFDFA1 function

scale = [16, 32, 64, 128, 256]; %Note: power of 2

% DFA looks at q=2, with MFDFA we can look at negative q-order giving more
% information about the structure.

q2= -2:0.1:2; % Second statistical moment
q3= -3:0.1:3; % Third statistical moment
q5 = -5:0.1:5;% Fifth statistical moment

% Polynomial order for detrending
m = 1;

% Flag for output plot (1 = on, 0 = off)
plot_flag = 1; 

% MFDFA starts here
figure(2);
[Hq,tq,hq,Dq,Fq] = MFDFA1(StrideIntervals, scale,q2 , m,1);

% MFDFA width: Larger width values can be interpreted as greater number of
% patterns in the time series
width = max(hq)- min(hq);

% What to look for?
%	  - Look at the multifractal spectrum width. The width will increase as
%       the q-order increases, but they are limitations.
%
%	  - To perform MFDFA with higher q-order you need a lot of data points.
%	    For short time series a q-order that ranges between -3 to 3 should be
%	    enough.

%% Section 2: Full Analysis

% Perform MFDFA analysis on a full dataset, which includes three overground
% walking conditions for three participants.
% Participants were instructed to synchronize their right step to various
% metronomes: Pink noise and White noise. Prior to those trials, participants
% performed a self-paced trial for baseline.

clear;
clc;
    
% Select folder where data is:
my_directory = uigetdir(matlabroot, 'select folder where data is');
S = dir(fullfile(my_directory, '*.mat')); % Returns a structure of the path specified in base path

% Select folder where you want to save the figures:
my_figures = uigetdir(matlabroot, 'select folder where you want to save your figure');

%Create variables for summary table
id_cond= [];
width=[];
Hurst_exp= [];
mass_exp= [];

for p = 1:numel(S) %The p loop goes through stride intervals folder
    fullpath = S(p).name;
    
    name= strsplit(fullpath, '_');
    ID= name{1,1};
    cond= name{1,2};
    ID_cond= [ID, cond]; % Combine the id name and condition

    disp(ID_cond) % display ID and condition to keep track of the progress of the loop

    Path_File= fullfile(S(p).folder,fullpath);
    data= load (Path_File); % Read matlab file
    
    % Create a time variable: select the first column of the stride intervals
    % variable (SI).
    time= data.SI(:,1);
    
    % Select the second column of the stride intervals
    % variable (SI) to create a 'StrideIntervals' variable.
    StrideIntervals= data.SI(:,2);
    
    %% Visual inspection before analysis:
    % Plot and Inspect your data:
    %      It's important that you look at your graph before analysis.
    %      Make sure that you don't see any abnormal large peaks or other
    %      abnormalities.
    
    %   What to look for?
    %      At this point of the pre-processing, your data should look
    %      noise-like time series (refer to slide 25-26 of the
    %      power-point for best practices.
    
    figure(1);
    plot(StrideIntervals, 'r-')
    title('Walking Stride intervals time series', 'FontSize', 16)
    xlabel('Time (s)');
    ylabel('Time Interval Amplitude (s)');
    disp('Press any key to continue!')
    
    %% Set MFDFA Parameters for Analysis:
    % MFDFA1 function parameters (Details can be found in the MFDFA1
    % function)
    
    scale = [16, 32, 64, 128, 256];
    
    % DFA looks at q=2, with MFDFA we can look at negative q-order giving more
    % information about the structure.
    
    q2= -2:0.1:2; % Second statistical moment
    q3= -3:0.1:3; % Third statistical moment
    q5 = -5:0.1:5;% Fifth statistical moment
    
    % Polynomial order for detrending
    m = 1;
    
    % Flag for output plot (1 = on, 0 = off)
    plot_flag = 1;
    
    % MFDFA starts here
    figure(2);
    [Hq,tq,hq,Dq,Fq] = MFDFA1(StrideIntervals, scale,q5 , m,1);
    
    % Export and save our figures into the folder titled "ANALYSIS OUTPUT"
    image= fullfile(my_figures, append(ID_cond, '.png'));
    saveas(gcf, image);
    
    disp('Press any key to continue!')
    
    % Hq: General Hurst exponent
    % tq: Mass exponent
    MFDFA_width= max(hq)- min(hq);
    Hq_avg= mean(Hq);
    tq_avg= mean(tq);
    
    
    %% Append all the variables in one table
    ID_cond= cellstr(ID_cond)
    id_cond= [id_cond;ID_cond,];
    width = [width; MFDFA_width,];
%     Hurst_exp= [Hurst_exp; Hq_avg,];
%     mass_exp= [mass_exp; tq_avg,];

end

%SUMMARY TABLE
T = table(id_cond,width);
% T = table(ID_cond, width, Hurst_exp, mass_exp);

%Save subject results
writetable(T,'MFDFA_ResultsSummary.csv')


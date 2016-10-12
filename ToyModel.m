clear all;
close all;
clc;

%%case with 1Ag
n_Ag = 1;
n_founders = 3;
rep = 9;
n_max_Bcells = n_founders*2^rep;
n_cycle_max = 200;
n_trial_max = 10;
a_act = 10;
a_threshold = 20;
a_min = -6;
conc = 0.9; %concentration doesn't induce frustration for conc [0.9, 1]

p_mut = 0.2; %per division.
p_CDR = 0.3;
p_FR_lethal = 0.9;
p_recycle = 0.85;
t_cell_selection = 0.6;

founders = rand(1,n_founders);
b_cells = zeros(n_trial_max, n_max_Bcells);
exit_cells = zeros(n_trial_max, n_cycle_max, floor(n_max_Bcells/4));
number_recycled_b_cells = zeros(n_trial_max, n_cycle_max);
number_exit_cells = zeros(n_trial_max, n_cycle_max);

tic;

%%INITIALIZATION + proliferation: 
%% each founder needs to meet a_act for the single Ag. b_cells all have affinity [a_act, a_act +1 ] at time 0.
cycle_number = 1;
number_recycled_b_cells(:,cycle_number) = n_founders;
number_exit_cells(:,cycle_number) = 0;

cycle_number = cycle_number +1;

for f = 1:n_founders
    f_start = (f-1)*2^rep+1;
    for b = f_start:f_start+2^rep-1
        b_cells(:,b) = founders(f) + a_act;
    end
end

number_recycled_b_cells(:,cycle_number) = size(b_cells,2);
number_exit_cells(:,cycle_number) = 0;

cycle_number = cycle_number +1;

%%%% The GC starts at cycle 3

%% Stochastic process: reproduce the GC reaction many times.

trial_number = 1;
while trial_number < n_trial_max
    disp(['TRIAL NUMBER ' num2str(trial_number)]);
    b_cells_trial = b_cells(trial_number,:);
    number_recycled_b_cells_trial = number_recycled_b_cells(trial_number,:);
    exit_cells_trial = exit_cells(trial_number, :, :);
    number_exit_cells_trial = number_exit_cells(trial_number,:);
    
    [b_cells_trial, number_recycled_b_cells_trial, exit_cells_trial, number_exit_cells_trial, final_cycle] = runTrial(b_cells_trial, exit_cells_trial, number_recycled_b_cells_trial, number_exit_cells_trial, conc, a_act, a_threshold, p_mut, p_CDR, p_FR_lethal, p_recycle, t_cell_selection, n_max_Bcells, n_cycle_max);
    
    for i = 1:size(b_cells_trial,2)
        b_cells(trial_number, i) = b_cells_trial(i);
    end   
    
    for i = 1: final_cycle
        number_recycled_b_cells(trial_number,i) = number_recycled_b_cells_trial(i);
        number_exit_cells(trial_number,i) = number_exit_cells_trial(i);
        
        for j = 1:size(exit_cells_trial,3)
            exit_cells(trial_number,i, j) = exit_cells_trial(1,i,j);
        end
    end

    trial_number = trial_number +1;
end

toc; 
%% Analyze trials
[pop_time, total_exit_cells, neutralized, breadth] = analysis( number_recycled_b_cells, number_exit_cells, exit_cells, n_trial_max, a_act, n_cycle_max, p_mut, p_recycle, t_cell_selection);
%function [ pop_time ] = analysis( number_recycled_b_cells, exit_cells, n_trial_max, a_act, n_cycle_max, p_mut, p_recycle, t_cell_selection )

%[breadth averageAffinity maxAffinity pop_over_time] = analysis(
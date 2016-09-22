function [ mutant ] = mutation( wildtype, a_act, a_threshold, p_CDR, p_FR_lethal)
%Takes a cell and mutates it. The type of mutation is determined and its 
%effect on the affinities for the two Ags is reported in the modified 
%mutant matrix.

p_CDR_lethal = 0.3;
p_CDR_silent = 0.5;
rand_CDR = rand;
rand_type = rand;
k = -0.7;
sigma = 1.3;
mu = -1.5;

mutant = wildtype;

%% CASE ALL MUTATIONS ARE IN CDR
rand_CDR = 0;

if rand_CDR < p_CDR
    %disp('case of mutation in CDR');
    if rand_type < p_CDR_lethal
        mutant = [];
    else if rand_type > p_CDR_lethal + p_CDR_silent
            mutant = wildtype + gevrnd(k, sigma, mu);
        end
    end
else
    disp('case of mutation in FR');
    if rand_type < p_FR_lethal
        %lethal FR mutation
        mutant = [];
    else
        if wildtype < a_act
            mutant = a_act;
        else if wildtype > a_threshold
                mutant = a_threshold;
            end
        end
    end    
end


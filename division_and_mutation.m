function [ daughters ] = division_and_mutation( b_cells_trial, a_act, a_threshold, p_mut, p_CDR, p_FR_lethal )

%division
daughters = [b_cells_trial b_cells_trial];
%round of SHM

for n = 1:size(daughters,2)
    rand_mut = rand;   
    if rand_mut < p_mut %if there is mutation, then call function mutation and modify the daughter cell
        %disp('Mutation ');
        %disp (['number of daughters before mutation ' num2str(size(daughters,2))]);
        mutant = mutation(daughters(n), a_act, a_threshold, p_CDR, p_FR_lethal);
        if ~isempty(mutant)
            daughters(n) = mutant;
        else            
            daughters(n) = NaN;
        end
    end
end
daughters = daughters(isnan(daughters) < 1);



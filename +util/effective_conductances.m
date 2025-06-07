% a Wrapper function to calculate effective conductances 
% allows to pass in a struct with multiple currents

function [g_eff_struct] = effective_conductances(all_results, vararg)
    % Initialize the output structure
    g_eff_struct = struct();
    for i = 1:length(vararg)
        name = vararg{i};
        E_ion = util.map_current_reversalPotential(name);
        current_field = ['I_' name];
        g_eff_struct.(name) = util.effective_conductance(all_results.(current_field), all_results.V, E_ion);
    end
end

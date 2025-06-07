function all_results = store_results(all_results, results, vararg)

 for i = 1:length(vararg)
     if isfield(results, vararg{i})
         all_results.(vararg{i}) = [all_results.(vararg{i}), results.(vararg{i})];
     else
         warning(['Field ', vararg{i}, ' not found in results.']);
     end
  end
end

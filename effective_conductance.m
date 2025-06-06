function g_ion = effective_conductance(I_ion, V, E_ion)
  g_ion = I_ion ./ (V - E_ion);
end

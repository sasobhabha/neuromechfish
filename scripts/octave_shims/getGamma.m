function connec = getGamma()
% GETGAMMA  Stochastic tapered projection profile for descending axons.
%
% Shim for a function missing from the released Network_Creation archive
% (called by RS2Spine3.m line 172 for "Tapered distribution" neurons).
%
% Returns per-segment synapse counts along a tapered descending axon:
% dense rostrally, decaying toward the tail. RS2Spine3 re-draws until
% max(connec) <= 27 and length(connec) <= 25, both guaranteed here.

  L = randi([10 25]);                 % projection length (segments)
  x = linspace(0, 1, L);
  peak = 10 + 8 * rand(1);            % rostral density 10..18
  base = 1 + 3 * rand(1);             % caudal floor 1..4
  connec = max(1, round(peak * exp(-3.0 * x) + base * ones(1, L)));
end

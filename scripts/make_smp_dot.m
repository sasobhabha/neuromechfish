% Generate the zebrafish Sensory-Motor Pathway (SMP) DOT graph using the
% original Peterson & Stobb MATLAB code, executed under GNU Octave.
% Output: _work/Network_Creation/SMP_generated.dot
%
% Mirrors makeSMP.m but skips drawDot.m (Graphviz/MATLAB plotting) --
% the DOT file is parsed in Python instead (scripts/build_connectome.py).

driver_dir = fileparts(mfilename('fullpath'));
smp_dir = fullfile(driver_dir, '..', '_work', 'Network_Creation');
addpath(smp_dir);
addpath(fullfile(driver_dir, 'octave_shims'));  % getGamma shim
cd(smp_dir);

% Seed RNGs for reproducible stochastic wiring
rand('state', 42); randn('state', 42);
more off;

% 1. Seed with the lateral-line portion of the graph
fid = fopen('lateral_line_stoch.dot', 'r');
txt = fscanf(fid, '%c');
fclose(fid);

out = 'SMP_generated.dot';
fid = fopen(out, 'w');
fprintf(fid, '%s', txt);
fclose(fid);
fprintf('1/5 lateral-line seed written\n'); fflush(stdout);

t0 = tic;

% 2. Brain: reticulospinal + brain interneurons
makeRSdot2(out);
fprintf('2/5 makeRSdot2 done (%.1fs)\n', toc(t0)); fflush(stdout);

% 3. Lateral line -> brain connections
LLtoRS(out, 1);
fprintf('3/5 LLtoRS done (%.1fs)\n', toc(t0)); fflush(stdout);

% 4. Descending neurons -> spinal cord
RS2Spine3(out, 1);
fprintf('4/5 RS2Spine3 done (%.1fs)\n', toc(t0)); fflush(stdout);

% 5. Spinal cord (interneurons + motor neurons)
makeSpinedot(out);
fprintf('5/5 makeSpinedot done (%.1fs)\n', toc(t0)); fflush(stdout);

% Close the digraph
fid = fopen(out, 'a');
fprintf(fid, '\n}\n');
fclose(fid);
fprintf('wrote %s\n', out); fflush(stdout);

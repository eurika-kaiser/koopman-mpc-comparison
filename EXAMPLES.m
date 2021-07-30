
clear all
close all

% Set path for saving figures
pathfigs = '../FIGURES/'; mkdir(pathfigs)
show_figs_publ = 1;

% Run model performance comparison for several systems
select_system = 'motor';
run_example(select_system, pathfigs, show_figs_publ)

select_system = 'duffing';
run_example(select_system, pathfigs, show_figs_publ)

select_system = 'van_der_pol';
run_example(select_system, pathfigs, show_figs_publ)

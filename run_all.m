% To do? This script assumes
% we are in the top level of this directory
main_paper = 1;
appendix = 1;

if main_paper
    fedput_calibration;
    optnomap_constpi_friedman;
    optnomap_comparison;
    optnomap_constpi_friedman;
    lc_icc_example_plots;
    OptMP_vary_alpha;
end

if appendix
    cd('appendix');
    fedput_calibration_constpi;
    optnomap_vary_lambda;
    optnomap_comparison_lambda;
    optnomap_vary_iput;
    optnomap_vary_ilaw;
    optnomap_vary_strike;
    optnomap_constrates;
    map_comparison_well_targeted;
    stability_fraction_03;
    stability_fraction_07;
    lc_poorly_targeted_example_plots;
    map_comparison_poorly_targeted;
    cd('bank_misalloc');
    fedput_calibration;
    optnomap_vs_no_bank_misalloc;
    cd('../../');
end

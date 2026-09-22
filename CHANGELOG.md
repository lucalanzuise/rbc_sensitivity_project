# Changelog

## Updated research version

This revision implements the methodological and reporting changes identified in
the project review:

- added a labour-disutility scale parameter (`chi`);
- fixed steady-state employment at `n_target = 1/3` across `sigma` and `varphi`
  experiments by recalibrating `chi`;
- made the `rho_z` and `rho_phi` grids symmetric at `{0.5, 0.7, 0.9, 0.98}`;
- added impact, peak, time-to-peak, half-life, cumulative-20Q and cumulative-40Q
  tables for every core sensitivity experiment;
- added a steady-state employment control table;
- added a direct productivity-versus-labour-disutility shock figure;
- added stacked variance-decomposition figures for the core and fiscal models;
- labelled theoretical moments explicitly as HP-filtered (`lambda = 1600`);
- exported autocorrelations through lag 5 where available;
- added a research-question-to-output map;
- added explicit balanced-budget lump-sum tax accounting to the fiscal model;
- documented the interpretation limits of calibrated variance decompositions.

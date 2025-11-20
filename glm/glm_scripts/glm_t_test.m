% GLM Arrays
p_pos_WT   = [0.0845, 1.17e-09, 7.51e-11, 0.0027, 0.249, 2.14e-09, 2.08e-06, 0.0222, 0, 6.47e-13, 0.0105, 0.206, 1.45e-07];
p_elem_WT  = [0.000224, 0, 0, 2.01e-12, 0.00305, 0, 0, 4.42e-09, 0, 0, 0.00252, 0.000197, 0];

R2_pos_WT  = [0.008, 0.054, 0.061, 0.018, 0.005, 0.053, 0.036, 0.012, 0.107, 0.072, 0.014, 0.006, 0.042];
R2_elem_WT = [0.024, 0.178, 0.099, 0.069, 0.017, 0.190, 0.214, 0.051, 0.309, 0.100, 0.018, 0.024, 0.129];

w_pos_WT   = [0.26, 0.23, 0.38, 0.20, 0.23, 0.22, 0.14, 0.19, 0.26, 0.42, 0.44, 0.19, 0.25];
w_elem_WT  = [0.74, 0.77, 0.62, 0.80, 0.77, 0.78, 0.86, 0.81, 0.74, 0.58, 0.56, 0.81, 0.75];

p_pos_KO   = [0.000286, 0.373, 0.14, 0.528, 0.132];
p_elem_KO  = [2.45e-09, 0.0141, 2.59e-08, 1.22e-11, 0.33];

R2_pos_KO  = [0.023, 0.004, 0.007, 0.003, 0.007];
R2_elem_KO = [0.052, 0.013, 0.047, 0.065, 0.004];

w_pos_KO   = [0.31, 0.23, 0.13, 0.04, 0.62];
w_elem_KO  = [0.69, 0.77, 0.87, 0.96, 0.38];

% Run t-tests for each measure
[h_wpos, p_wpos]     = ttest2(w_pos_WT,   w_pos_KO,   'Vartype','unequal');
[h_welem, p_welem]   = ttest2(w_elem_WT,  w_elem_KO,  'Vartype','unequal');
[h_r2pos, p_r2pos]   = ttest2(R2_pos_WT,  R2_pos_KO,  'Vartype','unequal');
[h_r2elem, p_r2elem] = ttest2(R2_elem_WT, R2_elem_KO, 'Vartype','unequal');

[h_ppos, p_ppos]     = ttest2(p_pos_WT,   p_pos_KO,   'Vartype','unequal');
[h_pelem, p_pelem]   = ttest2(p_elem_WT,  p_elem_KO,  'Vartype','unequal');

[h_r2wt, p_r2wt] = ttest2(R2_elem_WT, R2_pos_WT, 'Vartype','unequal');
[h_r2ko p_r2ko] = ttest2(R2_elem_KO, R2_pos_KO, 'Vartype','unequal');

fprintf('w_pos:   p=%.4f, significant=%d\n',  p_wpos,  h_wpos);
fprintf('w_elem:  p=%.4f, significant=%d\n',  p_welem, h_welem);
fprintf('R2_pos:  p=%.4f, significant=%d\n',  p_r2pos, h_r2pos);
fprintf('R2_elem: p=%.4f, significant=%d\n',  p_r2elem, h_r2elem);
fprintf('p_pos:   p=%.4f, significant=%d\n',  p_ppos,  h_ppos);
fprintf('p_elem:  p=%.4f, significant=%d\n',  p_pelem, h_pelem);


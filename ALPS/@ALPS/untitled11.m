close all
clear all
clc

f.filename='ED_ladder_2x10_Jl_1.34_Jr_3.42_H_Dep/ED_ladder_2x10_Jl_1.34_Jr_3.42_H_0.5';
f.tasks=2;
f.g=2.06;
f.measurements={'magnetization'};

o=ALPS_xml(f);

f1.filename='ED_ladder_2x10_Jl_1.34_Jr_3.42_T_Dep/ED_ladder_2x10_Jl_1.34_Jr_3.42_T_0.31';
f1.tasks=2;
f1.g=2.06;
f1.measurements={'magnetization','entropy'};

o1=ALPS(f1);
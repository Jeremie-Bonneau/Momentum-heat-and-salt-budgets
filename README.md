README for momentum, heat and salt budgets
Author: Jérémie Bonneau, UBC, September 22 2024
The files here are to reproduce the results by Bonneau et al. 2024. Journal of Physical Oceanography (submitted)
For any questions, please email me at jbonneau@mail.ubc.ca

What you need to do for the heat and salt budgets:
1. In salt_budget.m, heat_budget.m and coefficients.m, change the directories where the data (mis_hourly_avg.m, mel_hourly_avg.m, adcp.m) is loaded from to match your directory
2. Run salt_budget.m and save K
3. Run heat_budget.m and save Hmis and m
4. Run coefficients, you will need to point the function melt_3eq.m for this to run

What you need to do for the momentum budget
1. Change the directory where the ctd profiles are imported from (Milne2015Jul_ctdd.mat)
2. Run channel_ctd_vel.m

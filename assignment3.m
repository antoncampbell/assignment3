%% Assignment 3
% By: Anton Campbell, 100975168
%% Question 1 - Monte-Carlo for Linear Voltage
% The Assignment 1 was re-run without a bottle neck. The applied voltage
% was set to 0.1 V.
%
% The simulation was run:
electron_box_3modes_A3(2)
pause(2);
%%
% 
% *(a)*
%
% The electric field (as calculated above) is $5.000\times 10^{5} V/m$.
%
% *(b)*
%
% The force on each electron (as calculated above) is $8.011\times 10^{-14} N$.
%
% *(c)*
%
% The accelation of each electron (as calculated above) is $3.382\times 10^{17} m/s^2$.
%
% Figure 1 shows some particle trajectories for the linear voltage 
% gradient. The particles scatter after different distances traveled. The
% paths are curved towards the right as the acceleration is to the right.
%
% *(d)*
% 
% The current formula is:
%
% $$\bar{v_{x}} nqW$$
%
% where $\bar{v_{x}}$ is the average velocity in the x-direction, $n$ is 
% the electron concentration, $q$ is the charge of an electron, and $W$ is
% the length of the y-direction boundary.
%
% Figure 2 shows the current over time. The current starts low. 
% It increases as the electron accelerate but stop at around $0.012 A$ 
% as the current is limited by the collisions and mean free path.
%
% *(e)*
%
% Figure 3 shows the density map at the end of the simulation.
% Figure 4 shows the temperature map at the end of the simulation. 
% The density
% and temperature are both relatively uniform as the boundary is periodic 
% in the x-direction. Acceleration is only in the x-direction to the right.


%% Question 2 - Finite Difference Potential Solution
% The Assignment 2 was re-run for the bottle neck case. The applied voltage
% was set to 0.8 V.
%
% The simulation was run:
box_potential_2a_A3_plot(0.8)
pause(2);


%%
% 
% *(a)*
%
% Figure 5 shows the voltage map with equipotential lines. The
% equipotential lines have closer spacing at the 
% "bottle-neck". Thus, the greatest voltage drop is across the 
% "bottle-neck". The low conductivity near the "bottle-neck" equates to a 
% higher resistance. Since current should be conserved, voltage drop 
% is higher for higher resistance according to the equation $V=IR$.
%
% *(b)*
%
% Figure 6 shows a map of electric field vectors. The vectors are all
% pointing to the right. The magnitude of the electric field is much higher
% inside the boxes because the voltage change is much quicker as seen in 
% Figure 5. 
%

%% Question 3 - Combined Simulations
% The Assignment 2 was combined with Assignment 1 for the bottle neck case.
% The applied voltage was set to 0.8 V.
%
electron_box_3modes_A3(3)

%%
% 
% *(a)*
%
% Figure 1 shows some particle trajectories when the bottle neck is 
% present. The paths are curved towards the right as the acceleration is 
% to the right.
%
% *(b)*
% 
% Figure 3 shows the density map at the end of the simulation.
% The density is high to left of the boundaries. The acceleration is to the
% right so electrons often get stuck on the left.
%
% *(c)*
% 
% There are a few possible next steps to make the simulation more accurate. 
% The number of electrons could be increased for the Monte-Carlo simulation.
% The mesh density could be increased for the finite difference method of
% the voltage solution.


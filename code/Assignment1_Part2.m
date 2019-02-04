%% ELEC 4700 - Assignment 1 Part 2 - Andrew Paul 100996250 - Would like the option of a meeting
% The second section of this assignment introduces random scattering as
% well as the elctrons having different initial velocities. THe initial
% velocities are set randomly follwoing a Maxwell-Boltzmann distribution
% which follows a Gaussian distribution with a different standard deviation
% dictated by the specifics of a Maxwell-Boltzmann equation. The particles
% are given some probability of randomly scattering dictated by:
%
% $$P_{scat} = 1 - e^{-\frac{dt}{\tau_{mn}}}$$
%
% A histogram is displayed below he code for 20 electrons but a better
% distribution is shown in the figures below the code where 10000 electrons
% populated the silicon. The 2D plot of electron trajectories only accounts
% for 20 electrons for run time purposes and it is easier to see their path
% with less population. Each time the particles scatter their velocity is
% changed and dictated by a Maxwell-Boltzmann distribution, therefore their
% temperature is changed. The average temperature of the particles is taken
% after each time step and plotted. The plot below the code only shows the
% averaging of 20 electrons where the additional plots show an average temperature
% vs. time plot of 10000 electrons.

% list of constants
m0 = 9.11e-31;
mn = 0.26*m0;
kB = 1.38e-23;
T = 300;

%region limits
xlim = 200e-9;
ylim = 100e-9;

% thermal velocity
vth = sqrt(2*kB*T/mn);

%initialize the number of electrons
num_electrons = 20;

% defining array for electrons (x postion, y position, angle, velocity)
electron = zeros(num_electrons, 4);

% the previous position of the electron (previous x position, previous y
% position)
electron_prev = zeros(num_electrons, 2);

%spacial step should be smaller than 1/100 of region size
time_step = xlim/vth/100;
time_total = time_step*500;
%num_step = time_total/time_step;

% used to make each electron a different colour
electron_colour = hsv(num_electrons);

% counter used to check temperature is constant
count = 0;

% scattering probability
Pscat = 1-exp(-time_step/0.2e-12);

%set an initial random postion and a fixed velocity for each electron
for i=1:num_electrons
    for j=1:4
        if(j==1)
            electron(i,j) = xlim*rand();
        elseif(j==2)
            electron(i,j) = ylim*rand();
        elseif(j==3)
            electron(i,j) = 2*pi*rand();
        else
            electron(i,j) = randn()*vth;
        end
    end
end

figure(3)
hist(electron(:,4))
title('Velocity Distribution')

% define a temperature and time array for plotting
temperature= zeros(time_total/time_step,1);
time = zeros(time_total/time_step,1);

% counter for mean collision time
collision_count = 0;

running_time = 0;

% velocity array used to calculated mean free path
velocity = zeros(time_total/time_step,1);

% update each electrons positon for each time step
for k=0:time_step:time_total
    avg_temp = 0;
    avg_velocity = 0;
    for m=1:num_electrons
        % allows electrons to pass through to the other side of the region
        %in the x-direction
        if (electron(m,1) >= xlim)
            electron(m,1) = 0;
            electron_prev(m,1) = 0;
        elseif (electron(m,1) <= 0)
            electron(m,1) = xlim;
            electron_prev(m,1) = xlim;
        end
        % electrons are reflected at the same angle if they strike the limits
        % of the region in the y-driection
        if ((electron(m,2) >= ylim) || (electron(m,2) <= 0))
            electron(m,3) = pi - electron(m,3);
            electron(m,4) = -electron(m,4);
        end
        
        % see if the particle scatters or not
        if(Pscat > rand())
            % scatters at a random angle
            electron(m,3) = 2*pi*rand();
            % new velocity for scattering - gaussian with some
            % MAXWELL-BOLTZMAN standard deviation
            vx_new = randn()*vth;
            vy_new = randn()*vth;
            v_new = sqrt(vx_new^2+vy_new^2);
            electron(m,4) = v_new;
            collision_count =+ 1;
        end
        
        avg_temp = avg_temp + (electron(m,4)^2)*mn/(2*kB);
        avg_velocity = avg_velocity + electron(m,4);
        
        %plot the movement of each electron
        if(k~=0)
            figure(1)
            plot([electron_prev(m,1),electron(m,1)],[electron_prev(m,2),electron(m,2)],'color',electron_colour(m,:))
            axis([0 xlim 0 ylim]);
        end

    end

    title('Electron movement: random scattering')
    xlabel('x-axis position (m)')
    ylabel('y-axis position (m)')
    hold on
    pause(0.001)

   % set the previous postion of the electron to the current electron
   %postion for the next itteration
   electron_prev(:,1) = electron(:,1);
   electron_prev(:,2) = electron(:,2);
   
   % set the electron postion to an updated position
   electron(:,1) = electron(:,1) + cos(electron(:,3)).*electron(:,4).*time_step;
   electron(:,2) = electron(:,2) + sin(electron(:,3)).*electron(:,4).*time_step;
   
   count = count +1;
   temperature(count,1) = avg_temp/num_electrons;
   time(count,1) = k + time_step;
   velocity(count,1) = avg_velocity;
     
end

mean_collision = time_total/collision_count;
avg_vth = 0;
for n=1:500
    avg_vth =+ velocity(n,1);
end
avg_vth = avg_vth/size(velocity,1);

MFP = avg_vth*mean_collision;
   
figure(2)
plot(time,temperature)
axis([0 time_total, 0 1100])
title('Temperature of electrons over time')
xlabel('time (s)')
ylabel('Temperature (K)')

%% Data found when populating region with 10000 electrons
% The following histogram shows the inital velocity distribution of each
% electron following a Maxwell-Boltzmann distribution.
%%
% 
% <<Histogram_part2.jpg>>
%
%This clearly follows Maxwell-Boltzmann distribution as we increase the
%population of the region.
%The tempereature plot shown below is taken using the very populated region
%of 10000 electrons and it is clear that they electrons initially increase
%in temperature but then their average maintains somewhat constant.

%%
% 
% <<Temp_plot_part2.jpg>>
% 
%Finally, the mean time between collisions and mean free path was
%recalcualted for 10000 electrons using the same method as part 1. The mean
%time between colisions was found to be 5.3485 ps which is slightly larger
%than in part 1 but the mean time was given in part 1 and was likely
%estimated to be too small as the electrons did not possess random
%scattering properties. The mean collision time is still in the order of ps
%confirming that the model is operating properly. The mean free path was
%calculated to be 1.038 x $10^{-5}$ m. This is again larger than part 1 but
%is proportional to the mean collision time so it is expected.






% function cold_plasma
clear all
minX = 0;
maxX = 2;

% minX = 0.45;
% maxX = 0.55;

N = 200; 
vmin= -1; vmax = 1;
% vmin= -0.1; vmax = 0.1;

epsilon = 0.05;
delta = 0.002; %0.002, 0.05
omega_0 = 1; 


% dt = 0.001;
t_final = 5;
alpha = zeros(1,N);
for i = 1:N
    alpha(i) = (i-0.5)/N;
end

x_1 = zeros(1,N);
v_1 = zeros(1,N);


for i = 1:N
    x_1(i) = alpha(i) + epsilon*sin(2*pi*alpha(i));
    v_1(i) = 0;
end


num_interval = N;

    
% Get reference answer using Euler and dt = 0.001 
ref_x = x_1;
ref_v = v_1;
T_euler = t_final/0.0001;
for step = 1:T_euler
   x_tt = x_dd(ref_x,N,omega_0,delta);
   for i = 1:N
       ref_x(i) = ref_x(i) + 0.0001 * ref_v(i);
       ref_v(i) = ref_v(i) + 0.0001 * x_tt(i);
   end
end

clear x_tt

% Use RK4 with dt = 0.5^0 to 0.5^5
for power = 0:4
    for i = 1:N
        x_1(i) = alpha(i) + epsilon*sin(2*pi*alpha(i));
        v_1(i) = 0;
    end

    dt(power+1) = 0.5^power;
    Nstep(power+1) = t_final/dt(power+1); 
    for step = 1:Nstep
        x_tt = x_dd(x_1,N,omega_0,delta);
        for i = 1:N
            k1_x = v_1(i);
            k1_v = x_tt(i);

            k2_x = v_1(i) + 0.5*dt(power+1)*k1_x;
            k2_v = x_tt(i) + 0.5*dt(power+1) * k1_v;

            k3_x = v_1(i) + 0.5*dt(power+1)*k2_x;
            k3_v = x_tt(i) + 0.5*dt(power+1) * k2_v;

            k4_x = v_1(i) + dt(power+1)*k3_x;
            k4_v = x_tt(i) + dt(power+1) * k3_v;

            x_1(i) = x_1(i) + dt(power+1)/6 * (k1_x+2*k2_x+2*k3_x+k4_x);
            v_1(i) = v_1(i) + dt(power+1)/6 * (k1_v+2*k2_v+2*k3_v+k4_v);
        end
    end

    % get error 
    Error(power+1) = 0;
    Ev(power+1) = 0;
    for i = 1:N
        Error(power+1) = Error(power+1) + abs(x_1(i)-ref_x(i));
        Ev(power+1) = Ev(power+1) + abs(v_1(i)-ref_v(i));
    end
    % Error(power+1) = sqrt(Error(power+1));
    % Ev(power+1) = sqrt(Ev(power+1));
end

figure(1);
loglog(dt,Error);
hold on 
xlabel('dt');ylabel('Error');
y = dt.^4;
loglog(dt,y);
hold off
legend('Error', 'slope 4');


figure(2);
loglog(dt,Ev);
hold on 
xlabel('dt');ylabel('Error');
loglog(dt,y);
hold off
legend('Error', 'slope 4');



function Efield = x_dd(x,particle_sum,omega_0,delta)
    for i = 1:particle_sum
        x(i) = mod(x(i),1);
    end
    Efield = zeros(1,particle_sum);
    for i = 1:particle_sum
        pho_bar = 0;
        a = 0;
        kernel = 0;
        for j = 1:particle_sum
            % Efield(i) = Efield(i) - k(x(i),active(j),delta)* omega_0*(1/num_active);  
            kernel = kernel - k(x(i),x(j),delta)* omega_0*(1/particle_sum);
            pho_bar = pho_bar + (k(1,x(j),delta) - k(0,x(j),delta)) * omega_0 *(1/particle_sum);
            a = a + (g(1,x(j),delta) - g(0,x(j),delta)) * omega_0 *(1/particle_sum);            
        end
        Efield(i) = Efield(i) + kernel+  pho_bar *(x(i)-0.5) - a;
    end
end



function weight = k(x,y,delta)
    weight = 1/2*(x-y)/((x-y)^2+delta^2)^0.5;
end

function green = g(x,y,delta)
    green = -0.5*((x-y)^2+delta^2)^0.5;
end



    
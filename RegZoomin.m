% function cold_plasma
clear all
minX = 0.45;
maxX = 0.55;

% allN = 25;
allN = [25, 50,100,200]; 
vmin= -0.1; vmax = 0.1;
epsilon = 0.05;
delta = 0.002;
% pho_bar = 1;
omega_0 = 1;

dt = 0.01;
t_final = 5;
Nstep = t_final/dt; 

method = 1;

for index = 1:length(allN)
    N = allN(index);
    dx = 1/N; 

    % Initialize alpha, x and t
    alpha = zeros(1,N);
    for i = 1:N
        alpha(i) = (i-0.5)/N;
    end
    x = zeros(1,N);
    v = zeros(1,N);
    for i = 1:N
        x(i) = alpha(i) + epsilon*sin(2*pi*alpha(i));
        v(i) = 0;
    end

    str = sprintf('t = %0.5g',0);
    figure(1); plot(x,v); 
    xlabel('x'); ylabel('v'); title(str); axis([ minX maxX vmin vmax])

    if method == 1
        for step = 1:Nstep
            x_tt = x_dd(x,N,omega_0,delta);
            for i = 1:N
                x(i) = mod(x(i) + dt*v(i),1);
                v(i) = v(i) + dt* x_tt(i);
            end
          str = sprintf('t = %0.5g',step*dt);
          figure(1); plot(x,v); 
          xlabel('x'); ylabel('v'); title(str); axis([ minX maxX vmin vmax])
      
            if step == t_final/dt
              figure(2);subplot(2,2,index); plot(x,v,'-o','MarkerSize',1.25);axis([ minX maxX vmin vmax])
              
              hold on
              z = linspace(minX,maxX,15);
              y = zeros(length(z),1);
              plot(z,y,'--b')
              
              time = step*dt;
              title(sprintf('N = %.2f, dt = %.2f', N,dt));
            end
        end
    end

    if method == 2
        for step = 1:Nstep
            x_tt = x_dd(x,N,omega_0,delta);
            for i = 1:N
                k1_x = v(i);
                k1_v = x_tt(i);
                k2_x = v(i) + 0.5*dt*k1_x;
                k2_v = x_tt(i) + 0.5*dt * k1_v;
                k3_x = v(i) + 0.5*dt*k2_x;
                k3_v = x_tt(i) + 0.5*dt * k2_v;
    
                k4_x = v(i) + dt*k3_x;
                k4_v = x_tt(i) + dt * k3_v;
                
                x(i) = x(i) + dt/6 * (k1_x+2*k2_x+2*k3_x+k4_x);
                v(i) = v(i) + dt/6 * (k1_v+2*k2_v+2*k3_v+k4_v);

            end
          str = sprintf('t = %0.5g',step*dt);
          figure(1); plot(x,v); 
          xlabel('x'); ylabel('v'); title(str); axis([ minX maxX vmin vmax])
      
            if step == t_final/dt
              figure(2);subplot(2,2,index); plot(x,v,'-o','MarkerSize',1.25);axis([ minX maxX vmin vmax])
              
              hold on
              z = linspace(minX,maxX,15);
              y = zeros(length(z),1);
              plot(z,y,'--b')
              
              time = step*dt;
              title(sprintf('N = %.2f, dt = %.2f', N,dt));
            end
        end
    end


end



function Efield = x_dd(x,particle_sum,omega_0,delta)
    for i = 1:particle_sum
        x(i) = mod(x(i),1);
    end
    Efield = zeros(1,particle_sum);

    for i = 1:particle_sum
        for j = 1:particle_sum
            Efield(i) = Efield(i) + k(x(i),x(j),delta)* omega_0*(1/particle_sum);             
        end
    end
end


function weight = k(x,y,delta)
    c_delta = sqrt(1+4*delta^2);
    weight = -c_delta/2*(x-y)/((x-y)^2+delta^2)^0.5+x-y;
end


    
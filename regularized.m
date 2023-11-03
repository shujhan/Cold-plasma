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

% method 1 is Euler, 2 is RK4
method = 2;

dt = 0.01;
t_final = 20;
Nstep = t_final/dt; 


d1 = 0.05; % chord length of the interval

% Initialize alpha, x and t

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

x_2 = zeros(1,N);
v_2 = zeros(1,N);

x_3 = zeros(1,N);
v_3 = zeros(1,N);

x_4 = zeros(1,N);
v_4 = zeros(1,N);

x_5 = zeros(1,N);
v_5 = zeros(1,N);

for i = 1:N
    x_2(i) = x_1(i) + 1;
    v_2(i) = v_1(i);

    x_3(i) = x_1(i) - 1;
    v_3(i) = v_1(i);

    x_4(i) = x_1(i) + 2;
    v_4(i) = v_1(i);

    x_5(i) = x_1(i) - 2;
    v_5(i) = v_1(i);

    x_6(i) = x_1(i) + 3;
    v_6(i) = v_1(i);

end


num_interval = N;

str = sprintf('t = %0.5g',0);
figure(1); 
% plot(x,v)
plot(x_1,v_1,'r')
hold on
plot(x_2,v_2,'b')
hold on
plot(x_3,v_3,'g')
hold on
plot(x_4,v_4,'k')
hold on
plot(x_5,v_5,'m')
hold off


xlabel('x'); ylabel('v'); title(str); axis([ minX maxX vmin vmax])


  figure(2); subplot(3,2,1); plot(x_1,v_1,'-r');
  hold on
  plot(x_2,v_2,'-b');
  hold on
  plot(x_3,v_3,'g')
  hold on
  plot(x_4,v_4,'k')
  hold on
  plot(x_5,v_5,'m')
  hold off

  axis([ minX maxX vmin vmax])
  title('t = 0');
  part = 1;
    
  % Nstep = 1;
    for step = 1:Nstep
        x_tt = x_dd(x_1,N,omega_0,delta);
        if method == 1
            for i = 1:N
                x_1(i) = x_1(i) + dt * v_1(i);
                v_1(i) = v_1(i) + dt * x_tt(i);
            end

            for i = 1:N
                x_2(i) = x_1(i)+1;
                v_2(i) = v_1(i);
            
                x_3(i) = x_1(i) - 1;
                v_3(i) = v_1(i);
            
                x_4(i) = x_1(i) +2;
                v_4(i) = v_1(i);
            
                x_5(i) = x_1(i) -2;
                v_5(i) = v_1(i);

                x_6(i) = x_1(i) + 3;
                v_6(i) = v_1(i);

            end
        end


        if method ==2
            for i = 1:N
                k1_x = v_1(i);
                k1_v = x_tt(i);
                k2_x = v_1(i) + 0.5*dt*k1_x;
                k2_v = x_tt(i) + 0.5*dt * k1_v;
                k3_x = v_1(i) + 0.5*dt*k2_x;
                k3_v = x_tt(i) + 0.5*dt * k2_v;

                k4_x = v_1(i) + dt*k3_x;
                k4_v = x_tt(i) + dt * k3_v;

                x_1(i) = x_1(i) + dt/6 * (k1_x+2*k2_x+2*k3_x+k4_x);
                v_1(i) = v_1(i) + dt/6 * (k1_v+2*k2_v+2*k3_v+k4_v);
            end

            for i = 1:N
                x_2(i) = x_1(i)+1;
                v_2(i) = v_1(i);

                x_3(i) = x_1(i) - 1;
                v_3(i) = v_1(i);
            
                x_4(i) = x_1(i) +2;
                v_4(i) = v_1(i);
            
                x_5(i) = x_1(i) -2;
                v_5(i) = v_1(i);

                x_6(i) = x_1(i) + 3;
                v_6(i) = v_1(i);
            end
        end

      str = sprintf('t = %0.5g',step*dt);
      figure(1); 
      % plot(x,v); 

      plot(x_1,v_1,'-or','MarkerSize',2)
      hold on
      plot(x_2,v_2,'-ob','MarkerSize',2)
      hold on
      plot(x_3,v_3,'-og','MarkerSize',2)
      hold on
      plot(x_4,v_4,'-ok','MarkerSize',2)
      hold on
      plot(x_5,v_5,'-om','MarkerSize',2)
      hold on
      z = linspace(minX,maxX,15);
      y = zeros(length(z),1);
      plot(z,y,'--b')
      hold off


      xlabel('x'); ylabel('v'); title(str); axis([ minX maxX vmin vmax])

        % if step == 125 || step == 250 || step == 375 || step == 500 
        if step == 400 || step == 800 || step == 1200 || step == 1600 || step == 2000
        % if step == 500 
          figure(2);subplot(3,2,1+part); 
          % plot(x,v,'-o','MarkerSize',1.25);
          plot(x_1,v_1,'r')
          hold on
          plot(x_2,v_2,'b')
          hold on
          plot(x_3,v_3,'g')
          hold on
          plot(x_4,v_4,'k')
          hold on
          plot(x_5,v_5,'m')
          hold off
          axis([ minX maxX vmin vmax])
          time = step*dt;
          title(sprintf('t = %.2f', time));
          part = part+1;
        end
    end



function Efield = x_dd(x,particle_sum,omega_0,delta)
    for i = 1:particle_sum
        x(i) = mod(x(i),1);
    end
    Efield = zeros(1,particle_sum);
    for i = 1:particle_sum
        for j = 1:particle_sum
            Efield(i) = Efield(i) - k(x(i),x(j),delta)* omega_0*(1/particle_sum);             
        end
    end
end


function weight = k(x,y,delta)
    c_delta = (1+4*delta^2)^0.5;
    weight = -c_delta/2*(x-y)/((x-y)^2+delta^2)^0.5+x-y;
end




    
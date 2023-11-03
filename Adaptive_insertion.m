% function cold_plasma
clear all

minX = 0;
maxX = 2;

N = 200; 
vmin= -1; vmax = 1;

epsilon = 0.05;
delta = 0.002; %0.05,0.002
omega_0 = 1;

% method 1 is Euler, 2 is RK4
method = 1;

dt = 0.01;
t_final = 20;
Nstep = t_final/dt; 

d1 = 0.05; % chord length of the interval

% Initialize alpha, x and t

%particle_num is num particles in one period
particle_num = 2*N+1;

alpha = zeros(1,particle_num);
for i = 1:particle_num
    alpha(i) = (i-1)/(2*N);
end

x_1 = zeros(1,particle_num);
v_1 = zeros(1,particle_num);


for i = 1:particle_num
    x_1(i) = alpha(i) + epsilon*sin(2*pi*alpha(i));
    v_1(i) = 0;
end

    % only plot active points 
        num_active = (particle_num-1)/2;
        allActive = zeros(1,num_active);
        all_v = zeros(1,num_active);
        count = 1;
        for i = 1:particle_num
            if (rem(i,2) == 0)
                allActive(count) = x_1(i);
                all_v(count) = v_1(i);
                count = count +1;
            end
        end


        for i = 1:num_active
            x_2(i) = allActive(i)+1;
            v_2(i) = all_v(i);
            
            x_3(i) = allActive(i) - 1;
            v_3(i) = all_v(i);
            
            x_4(i) = allActive(i) +2;
            v_4(i) = all_v(i);
            
            x_5(i) = allActive(i) -2;
            v_5(i) = all_v(i);
        end


num_interval = N;

str = sprintf('t = %0.5g',0);
figure(1); 
% plot(x,v)
plot(allActive,all_v,'r')
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


  figure(2); subplot(3,2,1); 

  plot(allActive,all_v,'-r');
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
        % TODO: Only take active points into account to get x_tt
        x_tt = x_dd(x_1,particle_num,omega_0,delta);
        if method == 1
            for i = 1:particle_num
                x_1(i) = x_1(i) + dt * v_1(i);
                v_1(i) = v_1(i) + dt * x_tt(i);
            end
        end

        if method ==2
            for i = 1:particle_num
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
        end
    
    % only plot active points 
        num_active = (particle_num-1)/2;
        allActive = zeros(1,num_active);
        count = 1;
        for i = 1:particle_num
            if (rem(i,2) == 0)
                allActive(count) = x_1(i);
                all_v(count) = v_1(i);
                count = count +1;
            end
        end

        for i = 1:num_active
            x_2(i) = allActive(i)+1;
            v_2(i) = all_v(i);
            
            x_3(i) = allActive(i) - 1;
            v_3(i) = all_v(i);
            
            x_4(i) = allActive(i) +2;
            v_4(i) = all_v(i);
            
            x_5(i) = allActive(i) -2;
            v_5(i) = all_v(i);
        end

      str = sprintf('t = %0.5g',step*dt);
      figure(1); 
      % plot(x,v); 
      plot(allActive,all_v,'-r')
      hold on
      plot(x_2,v_2,'-b')
      hold on
      plot(x_3,v_3,'-g')
      hold on
      plot(x_4,v_4,'-k')
      hold on
      plot(x_5,v_5,'-m')
      hold on
      z = linspace(minX,maxX,15);
      y = zeros(length(z),1);
      plot(z,y,'--b')
      hold off

      xlabel('x'); ylabel('v'); title(str); axis([ minX maxX vmin vmax])

        % if step == 125 || step == 250 || step == 375 || step == 500 
        if step == 400 || step == 800 || step == 1200 || step == 1600 || step == 2000
          figure(2);subplot(3,2,1+part); 
          % plot(x,v,'-o','MarkerSize',1.25);
          plot(allActive,all_v,'r')
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




        % examine each interval in one period x_1, v_1
        count = 1;
        for i = 1:num_interval
            % euclidean distance in phase sapace
            dist = ((x_1(2*i+1) - x_1(2*i-1))^2 + (v_1(2*i+1) - v_1(2*i-1))^2)^0.5;
            if (dist > d1)
                % add the first point
                new_x_1(count) = x_1(2*i-1); 
                new_alpha(count) = alpha(2*i-1);
                new_v_1(count) = v_1(2*i-1);
                count = count +1;

                middle_alpha1 = (alpha(2*i-1) + alpha(2*i))/2;
                middle_alpha2 = (alpha(2*i+1) + alpha(2*i))/2;
                % quadratic interpolation using 3 points in the interval
                L01 = (middle_alpha1 - alpha(2*i))*(middle_alpha1 - alpha(2*i+1))/((alpha(2*i-1) - alpha(2*i))*(alpha(2*i-1) - alpha(2*i+1)));
                L11 = (middle_alpha1 - alpha(2*i-1))*(middle_alpha1 - alpha(2*i+1))/((alpha(2*i) - alpha(2*i-1))*(alpha(2*i) - alpha(2*i+1)));
                L21 = (middle_alpha1 - alpha(2*i-1))*(middle_alpha1 - alpha(2*i))/((alpha(2*i+1) - alpha(2*i-1))*(alpha(2*i+1) - alpha(2*i)));

                insert_x1 = x_1(2*i-1)*L01 + x_1(2*i)*L11+x_1(2*i+1)*L21;
                insert_v1 = v_1(2*i-1)*L01 + v_1(2*i)*L11+v_1(2*i+1)*L21;

                % add the second point (first insertion)
                new_x_1(count) = insert_x1;
                new_alpha(count) = middle_alpha1;
                new_v_1(count) = insert_v1;
                count = count +1;

                % add the third point 
                new_x_1(count) = x_1(2*i);
                new_alpha(count) = alpha(2*i);
                new_v_1(count) = v_1(2*i);
                count = count +1;
                L02 = (middle_alpha2 - alpha(2*i))*(middle_alpha2 - alpha(2*i+1))/((alpha(2*i-1) - alpha(2*i))*(alpha(2*i-1) - alpha(2*i+1)));
                L12 = (middle_alpha2 - alpha(2*i-1))*(middle_alpha2 - alpha(2*i+1))/((alpha(2*i) - alpha(2*i-1))*(alpha(2*i) - alpha(2*i+1)));
                L22 = (middle_alpha2 - alpha(2*i-1))*(middle_alpha2 - alpha(2*i))/((alpha(2*i+1) - alpha(2*i-1))*(alpha(2*i+1) - alpha(2*i)));

                insert_x2 = x_1(2*i-1)*L02 + x_1(2*i)*L12+x_1(2*i+1)*L22;
                insert_v2 = v_1(2*i-1)*L02 + v_1(2*i)*L12+v_1(2*i+1)*L22;

                % add the fourth point (second insertion)
                new_x_1(count) = insert_x2;
                new_alpha(count) = middle_alpha2;
                new_v_1(count) = insert_v2;
                count = count +1;     

            else 
                new_x_1(count) = x_1(2*i-1);
                new_alpha(count) = alpha(2*i-1);
                new_v_1(count) = v_1(2*i-1);
                count = count + 1;

                new_x_1(count) = x_1(2*i);
                new_alpha(count) = alpha(2*i);
                new_v_1(count) = v_1(2*i);
                count = count + 1;
            end
        end

        % add the last point
        new_x_1(count) = x_1(particle_num);
        new_alpha(count) = alpha(particle_num);
        new_v_1(count) = v_1(particle_num);

        % reset quatrature
        % count is the total particle number, which = 2*particle_num-1

        particle_num = count;
        num_interval = (particle_num-1)/2;

        x_1 = new_x_1;
        alpha = new_alpha;
        v_1 = new_v_1;
        clear new_x_1 new_v_1 new_alpha

        
    end



function Efield = x_dd(x,particle_sum,omega_0,delta)
    % for i = 1:particle_sum
    %     x(i) = mod(x(i),1);
    % end
    Efield = zeros(1,particle_sum);

    % only use active points with even index
    num_active = (particle_sum-1)/2;
    
    active = zeros(1,num_active);
    count = 1;
    for i = 1:particle_sum
        if (rem(i,2) == 0)
            active(count) = x(i);
            count = count +1;
        end
    end

    for i = 1:particle_sum
        for j = 1:num_active
            Efield(i) = Efield(i) + k(x(i),active(j),delta)* omega_0*(1/num_active);             
        end
    end
end


function weight = k(x,y,delta)
    c_delta = (1+4*delta^2)^0.5;
    weight = -c_delta/2*(x-y)/((x-y)^2+delta^2)^0.5+x-y;
end



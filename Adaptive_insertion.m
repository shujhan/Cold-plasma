% function cold_plasma
clear all

minX = 0;
maxX = 2;

N = 200; 
vmin= -1; vmax = 1;

epsilon = 0.05;
delta = 0.05; %0.05,0.002
omega_0 = 1;

% method 1 is Euler, 2 is RK4
method = 2;

dt = 0.04;
t_final = 20;
Nstep = t_final/dt; 

d1 = 0.05; % chord length of the interval

% passive points
alpha_passive = zeros(1,N+1);
x_passive = zeros(1,N+1);
v_passive = zeros(1,N+1);
% active points
alpha = zeros(1,N);
x = zeros(1,N);
v = zeros(1,N);
% weight
wt = zeros(1,N);

for i = 1:N+1
    alpha_passive(i) = (i-1)/N;
    x_passive(i) = alpha_passive(i) + epsilon*sin(2*pi*alpha_passive(i));
    v_passive(i) = 0;
end

for i = 1:N
    alpha(i) = 0.5*(alpha_passive(i)+alpha_passive(i+1));
    x(i) = alpha(i) + epsilon*sin(2*pi*alpha(i));
    v(i) = 0;
end

for i = 1:N
    wt(i) = alpha_passive(i+1) - alpha_passive(i);
end


% only plot active points 
x_2 = zeros(1,N);
v_2 = zeros(1,N);
% x_3 = zeros(1,N);
% v_3 = zeros(1,N);
% x_4 = zeros(1,N);
% v_4 = zeros(1,N);
% x_5 = zeros(1,N);
% v_5 = zeros(1,N);
% x_6 = zeros(1,N);
% v_6 = zeros(1,N);
for i = 1:N
    x_2(i) = x(i)+1;
    v_2(i) = v(i);
            
    % x_3(i) = x(i)-1;
    % v_3(i) = v(i);
    % 
    % x_4(i) = x(i)+2;
    % v_4(i) = v(i);

    % x_5(i) = x(i)-2;
    % v_5(i) = v(i);
    % 
    % x_6(i) = x(i)+3;
    % v_6(i) = v(i);
end

num_interval = N;

str = sprintf('t = %0.5g',0);
figure(1); 
plot(x,v,'r')
hold on
plot(x_2,v_2,'b')
hold on
% plot(x_3,v_3,'g')
% hold on
% plot(x_4,v_4,'k')
% hold on
% plot(x_5,v_5,'m')
% hold on
% plot(x_6,v_6,'m')



xlabel('x'); ylabel('v'); title(str); axis([ minX maxX vmin vmax])


  figure(2); subplot(3,2,1); 

  plot(x,v,'-r');
  hold on
  plot(x_2,v_2,'-b');
  hold on
  % plot(x_3,v_3,'g')
  % hold on
  % plot(x_4,v_4,'k')
  % hold on
  % plot(x_5,v_5,'m')
  % hold on
  % plot(x_6,v_6,'m')
  % hold off

  axis([ minX maxX vmin vmax])
  title('t = 0');
  part = 1;
    


    for step = 1:Nstep
        % TODO: Only take active points into account to get E field
        active_num = length(x);
        passive_num = length(x_passive);

        E_field = x_dd(x,x,omega_0,delta,wt);
        E_field_passive = x_dd(x_passive(1:passive_num-1),x,omega_0,delta,wt);
        

        if method == 1
            for i = 1:active_num
                x(i) = x(i) + dt * v(i);
                v(i) = v(i) + dt * E_field(i);
                x_passive(i) = x_passive(i) + dt * v_passive(i);
                v_passive(i) = v_passive(i) + dt * E_field_passive(i);
            end
            x_passive(passive_num) = x_passive(1) + 1;
            v_passive(passive_num) = v_passive(1);
        end

        if method ==2
            for i = 1:active_num
                k1_x = v(i);
                k1_v = E_field(i);
                k2_x = v(i) + 0.5*dt*k1_x;
                k2_v = E_field(i) + 0.5*dt * k1_v;
                k3_x = v(i) + 0.5*dt*k2_x;
                k3_v = E_field(i) + 0.5*dt * k2_v;
                k4_x = v(i) + dt*k3_x;
                k4_v = E_field(i) + dt * k3_v;
                x(i) = x(i) + dt/6 * (k1_x+2*k2_x+2*k3_x+k4_x);
                v(i) = v(i) + dt/6 * (k1_v+2*k2_v+2*k3_v+k4_v);
            end

            for i = 1:active_num
                k1_x = v_passive(i);
                k1_v = E_field_passive(i);
                k2_x = v_passive(i) + 0.5*dt*k1_x;
                k2_v = E_field_passive(i) + 0.5*dt * k1_v;
                k3_x = v_passive(i) + 0.5*dt*k2_x;
                k3_v = E_field_passive(i) + 0.5*dt * k2_v;
                k4_x = v_passive(i) + dt*k3_x;
                k4_v = E_field_passive(i) + dt * k3_v;
                x_passive(i) = x_passive(i) + dt/6 * (k1_x+2*k2_x+2*k3_x+k4_x);
                v_passive(i) = v_passive(i) + dt/6 * (k1_v+2*k2_v+2*k3_v+k4_v);
            end
            x_passive(passive_num) = x_passive(1) + 1;
            v_passive(passive_num) = v_passive(1);
        end
    
    % only plot active points 
        for i = 1:active_num
            x_2(i) = x(i)+1;
            v_2(i) = v(i);
            
            % x_3(i) = x(i) - 1;
            % v_3(i) = v(i);
            % 
            % x_4(i) = x(i) +2;
            % v_4(i) = v(i);
            
            % x_5(i) = x(i) -2;
            % v_5(i) = v(i);
            % 
            % x_6(i) = x(i) +3;
            % v_6(i) = v(i);
        end

      str = sprintf('t = %0.5g',step*dt);
      figure(1); 
      % plot(x,v); 
      plot(x,v,'-r')
      hold on
      % plot(x_2,v_2,'-b')
      % hold on
      % plot(x_3,v_3,'-g')
      % hold on
      % plot(x_4,v_4,'-k')
      % hold on
      % plot(x_5,v_5,'-m')
      % hold on
      % plot(x_6,v_6,'-m')
      % hold on
      z = linspace(minX,maxX,15);
      y = zeros(length(z),1);
      plot(z,y,'--b')
      hold off

      xlabel('x'); ylabel('v'); title(str); axis([ minX maxX vmin vmax])

        % if step == 125 || step == 250 || step == 375 || step == 500 
        if step == 400 || step == 800 || step == 1200 || step == 1600 || step == 2000
          figure(2);subplot(3,2,1+part); 
          plot(x,v,'r')
          hold on
          % plot(x_2,v_2,'b')
          % hold on
          % plot(x_3,v_3,'g')
          % hold on
          % plot(x_4,v_4,'k')
          % hold on
          % plot(x_5,v_5,'m')
          % hold on
          % plot(x_6,v_6,'-m')
          hold off
          axis([ minX maxX vmin vmax])
          time = step*dt;
          title(sprintf('t = %.2f', time));
          part = part+1;
        end




        % examine each interval in one period x,v
        % active_num equals interval numbers 
        count = 1;
        count_passive = 1;
        for i = 1:num_interval
            a0 = alpha_passive(i);
            a1 = alpha(i);
            a2 = alpha_passive(i+1);

            x0 = x_passive(i);
            x1 = x(i);
            x2 = x_passive(i+1);

            v0 = v_passive(i);
            v1 = v(i);
            v2 = v_passive(i+1);

            % euclidean distance in phase sapace
            dist = sqrt((x2-x0)^2 + (v2-v0)^2);
           
            if (dist > d1)
                % add first point to passive
                new_x_passive(count_passive) = x0;
                new_alpha_passive(count_passive) = a0;
                new_v_passive(count_passive) = v0;
                count_passive = count_passive + 1;
                
                % add second point to passive (active point on this interval)
                new_x_passive(count_passive) = x1;
                new_alpha_passive(count_passive) = a1;
                new_v_passive(count_passive) = v1;
                count_passive = count_passive + 1;


                % insert two active points new_alpha 
                middle_alpha1 = 0.5*(a1 + a0);
                middle_alpha2 = 0.5*(a1 + a2);
                % use newtons form interpolation
                a = ((x2-x1)/(a2-a1) - (x1-x0)/(a1-a0))/(a2-a0);
                b = (x1-x0)/(a1-a0);
                c = x0;
                insert_x1 = a*(middle_alpha1 - a0)*(middle_alpha1 - a1) + b*(middle_alpha1-a0) + c;
                insert_x2 = a*(middle_alpha2 - a0)*(middle_alpha2 - a1) + b*(middle_alpha2-a0) + c;
                
                a_v = ((v2-v1)/(a2-a1) - (v1-v0)/(a1-a0))/(a2-a0);
                b_v = (v1-v0)/(a1-a0);
                c_v = v0;
                insert_v1 = a_v*(middle_alpha1 - a0)*(middle_alpha1 - a1) + b_v*(middle_alpha1-a0) + c_v;
                insert_v2 = a_v*(middle_alpha2 - a0)*(middle_alpha2 - a1) + b_v*(middle_alpha2-a0) + c_v;
                
                new_alpha(count) = middle_alpha1;
                new_x(count) = insert_x1;
                new_v(count) = insert_v1;
                count = count + 1;

                new_alpha(count) = middle_alpha2;
                new_x(count) = insert_x2;
                new_v(count) = insert_v2;
                count = count + 1;
  

            else 
                new_x_passive(count_passive) = x_passive(i);
                new_alpha_passive(count_passive) = alpha_passive(i);
                new_v_passive(count_passive) = v_passive(i);
                count_passive = count_passive + 1;


                new_x(count) = x(i);
                new_alpha(count) = alpha(i);
                new_v(count) = v(i);
                count = count + 1;
            end
        end

        % add the last point
        new_x_passive(count_passive) = new_x_passive(1) + 1;
        new_alpha_passive(count_passive) = new_alpha_passive(1) + 1;
        new_v_passive(count_passive) = new_v_passive(1);

        % reset quatrature
        num_interval = length(new_alpha);

        x = new_x;
        alpha = new_alpha;
        v = new_v;

        x_passive = new_x_passive;
        alpha_passive = new_alpha_passive;
        v_passive = new_v_passive;

        clear new_x new_v new_alpha new_x_passive new_v_passive new_alpha_passive

        for i = 1:num_interval
            wt(i) = alpha_passive(i+1) - alpha_passive(i);
        end
    end



function Efield = x_dd(x_target,x_source,omega_0,delta,wt)
    target_num = length(x_target);
    source_num = length(x_source);
    for i = 1:source_num
        x_source(i) = mod(x_source(i),1);
    end
    for i = 1:target_num
        x_target(i) = mod(x_target(i),1);
    end


    Efield = zeros(1,target_num);

    for i = 1:target_num
        for j = 1:source_num
            Efield(i) = Efield(i) + k(x_target(i),x_source(j),delta)* omega_0*(wt(j));             
        end
    end
end


function kernel = k(x,y,delta)
    c_delta = (1+4*delta^2)^0.5;
    kernel = -c_delta/2*(x-y)/sqrt((x-y)^2+delta^2)+x-y;
end


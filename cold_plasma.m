% function cold_plasma
clear all
minX = 0;
maxX = 2;

N = 25; 
vmin= -1; vmax = 1;
epsilon = 0.05;
pho_bar = 1;
omega_0 = 1;

dt = 0.01;
t_final = 5;
Nstep = t_final/dt; 

dx = 1/N; 

% Initialize alpha, x and t
alpha = zeros(1,N);
for i = 1:N
    alpha(i) = (i-0.5)/N;
end
x = zeros(1,2*N);
v = zeros(1,2*N);
for i = 1:N
    x(i) = alpha(i) + epsilon*sin(2*pi*alpha(i));
    x(i+N) = x(i)+1;
    v(i) = 0;
    v(i+N) = v(i);
end

str = sprintf('t = %0.5g',0);
figure(1); plot(x,v); 
xlabel('x'); ylabel('v'); title(str); axis([ minX maxX vmin vmax])

  figure(2); subplot(3,2,1); plot(x,v,'-o','MarkerSize',1.25);axis([ minX maxX vmin vmax])
  title('t = 0');
  part = 1;



for step = 1:Nstep
    % for i = 1:N
    %     x_dd(i) = 0;
    %     a(i) = 0;
    %     for j = 1:N
    %         x_dd(i) = x_dd(i) - k(x(i),x(j))* omega_0*(1/N);
    %         a(i) = a(i) + x(j)*omega_0*(1/N);
    %     end
    %     a(i) = a(i)-1/2;
    %     x_dd(i) = x_dd(i) + pho_bar *(x(i)-0.5) - a(i);
    % end
    x_tt = x_dd(x,N,omega_0,pho_bar);
    for i = 1:N
        x(i) = mod(x(i) + dt*v(i),1);
        v(i) = v(i) + dt * x_tt(i);
        x(i+N) = mod(x(i+N) + dt*v(i+N),1)+1;
        v(i+N) = v(i+ N) + dt* x_tt(i);
    end
  str = sprintf('t = %0.5g',step*dt);
  figure(1); plot(x,v); 
  xlabel('x'); ylabel('v'); title(str); axis([ minX maxX vmin vmax])
  
    if step == 100 || step == 200 || step == 300 || step == 400 || step == 500
      figure(2);subplot(3,2,1+part); plot(x,v,'-o','MarkerSize',1.25);axis([ minX maxX vmin vmax])
      
      hold on
      z = linspace(minX,maxX,15);
      y = zeros(length(z),1);
      plot(z,y,'--b')
      time = step*dt;
      title(sprintf('t = %.2f', time));
      part = part+1;
    end
end




function Efield = x_dd(x,N,omega_0,pho_bar)
    % Efield = zeros(N,1);
    for i = 1:N
        Efield(i) = 0;
        a(i) = 0;
        for j = 1:N
            Efield(i) = Efield(i) - k(x(i),x(j))* omega_0*(1/N);
            a(i) = a(i) + x(j)*omega_0*(1/N);
        end
        a(i) = a(i)-1/2;
        Efield(i) = Efield(i) + pho_bar *(x(i)-0.5) - a(i);
    end
end


function weight = k(x,y)
  weight = 0;
  if x>y
      weight = 0.5;
  end
  if x < y
      weight = -0.5;
  end
end

    
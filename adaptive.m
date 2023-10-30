% function cold_plasma
clear all
minX = 0;
maxX = 2;


N = 200; 
vmin= -1; vmax = 1;
epsilon = 0.05;
delta = 0.04;
% pho_bar = 1;
omega_0 = 1;

tol = 0.05;



dt = 0.04;
t_final = 20;
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
    for i = 1:N
        x_dd = 0;
        a = 0; pho_bar = 0;
        for j = 2:2:(N-1)
            x_dd = x_dd - k(x(i),x(j),delta)* omega_0*(1/N);
            a = a+ ((-0.5*((1-x(j))^2+delta^2)^0.5)-(-0.5*((0-x(j))^2+delta^2)^0.5))*omega_0*(1/N);
            pho_bar = pho_bar + (0.5*(1-x(j))/((1-x(j))^2+delta^2)^0.5- 0.5*(0-x(j))/((0-x(j))^2+delta^2)^0.5)*omega_0*(1/N);
        end
        x_dd = x_dd + pho_bar *(x(i)-0.5) - a;
        
%       use rk4 with dt = 0.04

        v(i) = v(i) + dt* x_dd;
        v(i+N) = v(i+ N) + dt*x_dd;
        x(i) = x(i) + dt*v(i);
        x(i+N) = x(i+N) + dt*v(i+N);
    end
  str = sprintf('t = %0.5g',step*dt);
  figure(1); plot(x,v); 
  xlabel('x'); ylabel('v'); title(str); axis([ minX maxX vmin vmax])
  
%     if step == 100 || step == 200 || step == 300 || step == 400 || step == 500
    if step == 400 || step == 800 || step == 1200 || step == 1600 || step == 2000
      figure(2);subplot(3,2,1+part); plot(x,v,'-o','MarkerSize',1.25);axis([ minX maxX vmin vmax])
      time = step*dt;
      title(sprintf('t = %.2f', time));
      part = part+1;
    end
    
%   check if need insert
    for i = 1:2:N-2
        d1 = ((x(i)-x(i+2))^2+(v(i)-v(i+2))^2)^0.5;
        d2 = dist2(x(i),v(i),x(i+2),v(i+2),x(i+1),v(i+1));
        if d1 > tol || d2 > tol 
            N = N+2;
            p = polyfit(x(i:i+2),v(i:i+2),2);
 
            x_new = zeros(1,2*(N+2));
            x_new(1:i) = x(1:i);
            x_new(i+1) = abs(x(i+1)+x(i))/2;
            x_new(i+2) = x(i+1);
            x_new(i+3) = abs(x(i+1)+x(i+2))/2;
            x_new(i+4:N+2) = x(i+2:N);
            x_new(N+3:2*(N+2)) = x_new(1:N+2)+1;
            
            x = x_new;
            
            v_new = zeros(1,2*(N+2));
            v_new(1:i) = v(1:i);
            v_new(i+1) = polyval(p,abs(x(i+1)+x(i))/2);
            v_new(i+2) = v(i+1);
            v_new(i+3) = polyval(p,abs(x(i+1)+x(i+2))/2);
            v_new(i+4:N+2) = v(i+2:N);
            v_new(N+3:2*(N+2)) = v_new(1:N+2);
            
            v = v_new;
        end  
    end
end


function weight = k(x,y,delta)
    weight = 1/2*(x-y)/((x-y)^2+delta^2)^0.5;
end

function dist = dist2(x1,v1,x2,v2,x3,v3)
    numerator = abs((x2-x1)*(v1-v3)-(x1-x3)*(v2-v1));
    denominator = sqrt((x2-x1)^2+(v2-v1)^2);
    dist = numerator ./ denominator;
end

    
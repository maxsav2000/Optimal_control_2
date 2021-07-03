%% Задача нелинейного управления первая задача %Блок счета

%x=(v,m)
%\cdot x_1 = -g+u(v+l)/m
%\cdot x_2 = -u
t0 = 0; % начальное время
T = 2; % Конечное время

v0 = 0;
g = 1;
umax = 5;
m0 = 2;
M = 1;
x0=[v0, m0]';
eps = 0.1;
l = 2;


t_pers = [];
size_set = 200;
t_set = linspace(t0,T,size_set);
% стандартизированная сетка результатов

t_per1_vec = linspace(t0,T,2*size_set); % время для первого переключения
t_per2_vec = linspace(t0,T,size_set); % время для второго переключения

best_t  = [];
best_y = [];
%случай с одним переключением
ind_per = 1;

H_max = -1;
for i=2:2*size_set-1
    opts1 = odeset('RelTol',1e-5,'AbsTol',1e-5,'Refine',5,'MaxStep',1e-2);

    t_per = t_per1_vec(i);
    tspan1 = [t0 t_per];
    [t1,y1] = ode45(@(t,y) odefcn_1(t,y,umax,g,l), tspan1, x0, opts1);
    
    
    tspan2 = [t_per T];
    q = length(y1);
    x1 = [y1(q,1) y1(q,2)];
    [t2,y2] = ode45(@(t,y) odefcn_3(t,y,umax,g,l), tspan2, x1, opts1);
    
    q = length(y2);
    t = [t1' t2'];
    y = [y1' y2'];
    curr_H = trapz(t, y(1,:));
    if (abs(y2(q,1))<eps && y2(q,2)>M-0.01  && curr_H > H_max)
        H_max = curr_H;
        best_t  = t;
        best_y = y;
        abs(y2(q,1));
        ind_per = 1;
        t_pers = [t_per];
    end
end

for i=2:size_set-1
    t_per_1 = t_per2_vec(i);
    for j=i+1:size_set-1
        t_per_2= t_per2_vec(j);
        opts1 = odeset('RelTol',1e-3,'AbsTol',1e-3,'Refine',2,'MaxStep',1e-2);

        tspan1 = [t0 t_per_1];
        [t1,y1] = ode45(@(t,y) odefcn_1(t,y,umax,g,l), tspan1, x0, opts1);


        tspan2 = [t_per_1 t_per_2];
        q = length(y1);
        x1 = [y1(q,1) y1(q,2)];
        [t2,y2] = ode45(@(t,y) odefcn_2(t,y,umax,g,l,y1(q,1)), tspan2, x1, opts1);
        
        tspan3 = [t_per_2 T];
        q = length(y2);
        x2 = [y2(q,1) y2(q,2)];
        [t3,y3] = ode45(@(t,y) odefcn_3(t,y,umax,g,l), tspan3, x2, opts1);

        q = length(y3);
        t = [t1' t2' t3'];
        y = [y1' y2' y3'];
        curr_H = trapz(t, y(1,:));
        if (abs(y3(q,1))<eps && y3(q,2)>M-0.01 && curr_H > H_max)
            H_max = curr_H;
            best_t  = t;
            best_y = y;
            abs(y3(q,1));
            ind_per = 2;
            t_pers = [t_per_1 , t_per_2];
        end
    end
end%Случай с особым режимом


disp('Количество переключений:')
disp(ind_per);
disp('Моменты переключений:')
disp(t_pers);
disp('Максимальная высота:')
disp(H_max);

%% Графики m(t) v(t)

subplot(2,1,1);
plot(best_t, best_y(1,:),'b',[0,T],[eps,eps],'r',[0,T],[-eps,-eps],'r');
legend('v(t)','|v(t)|=eps');

xlabel('t');
ylabel('v(t)');
subplot(2,1,2);
plot(best_t, best_y(2,:),'b',[0,T], [M,M],'r',[0,T],[m0,m0],'k');
legend('m(t)','m=M','m=m_0');
xlabel('t');
ylabel('m(t)');
%Блок счета
%% Графики u(t) H(t)
subplot(2,1,1);
H_best = cumtrapz(best_t,best_y(1,:));
plot(best_t, H_best,'b');
legend('H(t)');

xlabel('t');
ylabel('H(t)');
subplot(2,1,2);

plot([0, t_pers(1),t_pers(1)+0.001,T], [umax umax,0,0]);
legend('u(t)');
xlabel('t');
ylabel('u(t)');




function dydt = odefcn_1(t,y,u_max,g,l)
dydt = zeros(2,1);
dydt(1) =-g+u_max*(y(1)+l)/y(2);
dydt(2) = -u_max;
end
function dydt = odefcn_2(t,y,u_max,g,l,v_c)
dydt = zeros(2,1);
dydt(1) = 0;
dydt(2) = -g*y(2)/(v_c+l);
end
function dydt = odefcn_3(t,y,u_max,g,l)
dydt = zeros(2,1);
dydt(1) = -g;
dydt(2) = 0;
end
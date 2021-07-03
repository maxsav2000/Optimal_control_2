%% Задача нелинейного управления вторая задача%Блок счета

t0 = 0;
v0 = 0;
T = 2;
M = 1;
m0 = 3;
umax = 1.4;
l = 2;
g = 1;
H = 0.7;
H_eps = 0.01;
alpha = 1;


size_set = 1000;
t_set = linspace(t0,T,size_set);
% стандартизированная сетка результатов
best_y = [];

J_min = umax^4 * T + alpha*umax*T+1; %минимизация функционала%
psi_N = 20;    
psi00 = linspace(-1, -1/psi_N, psi_N);
psi10 = linspace(-1, 1, 2 * psi_N);
psi30 = linspace(-1, 1, 2 * psi_N);

best_psi0 = 0;
best_psi1 = [];
best_psi2 = [];
best_psi3 = 0;
best_u = [];
best_v = [];
best_m = [];
best_t = [];

for i=1:(2 * psi_N)
    for j = 1:(2 * psi_N)
        for k = 1:(psi_N - 1)
            if (psi10(i)^2 + psi30(j)^2 + psi00(k)^2 > 1)
                continue;
            end
            psi20 = -sqrt(1 - psi10(i)^2 - psi30(j)^2 - psi00(k)^2);
            
            %\dot v = -g u (v+l)/m
            %\dot m = -u
            %\dot psi_1 = -psi_1 u/m - psi 3
            %\dot psi_2 = psi_1 u(v+l)/m^2
            
            
            % psi_0 = const, psi_3 = const
            % y = [v m psi_1 psi_2]
            tspan = [t0 T];
            x0 = [v0 m0 psi10(i) psi20];
            opts1 = odeset('RelTol',1e-5,'AbsTol',1e-5,'Refine',5,'MaxStep',1e-2);
            [t1,y1] = ode45(@(t,y) odefcn_1(t,y,umax,g,l,alpha,psi00(k),psi30(j),M), tspan, x0, opts1);
            h = trapz(t1,y1(:,1));
            
            if abs(H-h)<H_eps
                curr_u = zeros(1,length(t1));
                for q = 1:length(t1)
                    G = (-alpha*psi00(k) - y1(q,3)*(y1(q,1)+l)/y1(q,2)+y1(q,4))/(4*psi00(k));
                    u=0;
                    if (G > umax^3)
                        u = umax;
                    else
                        if (G < 0)
                            u = 0;
                        else
                            u = G^(1/3);
                        end
                    end
                    if (y1(q,2)<=M)
                        u = 0;
                    end
                    curr_u(q) = u;
                end
                f_u = curr_u.^4+alpha.*curr_u;
                curr_J = trapz(t1,f_u);
                if (curr_J<J_min)
                    J_min = curr_J;
                    best_t = t1;
                    best_psi0 = psi00(k);
                    best_psi1 = y1(:,3);
                    best_psi2 = y1(:,4);
                    best_psi3 = psi30(j);
                    best_u = curr_u;
                    best_v = y1(:,1);
                    best_m = y1(:,2);
                end
            end
        end
    end
end

if (J_min > umax^4 * T + alpha*umax)
    disp('Задача не разрешима:');
else
    disp('Минимальное значение функционала:');
    disp(J_min);
    disp('Начальные значения psi:');
    disp('psi00')
    disp(best_psi0);
    disp('psi10')
    disp(best_psi1(1));

    disp('psi20')
    disp(best_psi2(1));

    disp('psi30')
    disp(best_psi3);
end

%% Графики m(t) v(t)

subplot(2,1,1);
plot(best_t, best_v);
legend('v(t)');
xlabel('t');
ylabel('v(t)');
subplot(2,1,2);
plot(best_t, best_m,'b',[0,T], [M,M],'r',[0,T],[m0,m0],'k');
legend('m(t)','m=M','m=m_0');
xlabel('t');
ylabel('m(t)');
%% Графики u(t) H(t)
subplot(2,1,1);
H_best = cumtrapz(best_t,best_v);
plot(best_t, H_best,'b',[0,T],[H+H_eps,H+H_eps],'r',[0,T],[H-H_eps,H-H_eps],'r');
legend('H(t)','|H(t)-H|=H_{eps}');

xlabel('t');
ylabel('H(t)');
subplot(2,1,2);

plot(best_t, best_u,'b',[0,T],[umax,umax],'r');
legend('u(t)','u_{max}');
xlabel('t');
ylabel('u(t)');
%% Графики psi_1(t) psi2(t)
subplot(2,1,1);
plot(best_t, best_psi1,'b');
legend('psi_1(t)');

xlabel('t');
ylabel('psi_1(t)');
subplot(2,1,2);

plot(best_t, best_psi2,'b');
legend('psi_2(t)');
xlabel('t');
ylabel('psi_2(t)');





function dydt = odefcn_1(t,y,u_max,g,l,alpha,psi00, psi30, M)
dydt = zeros(4,1);
G = (-alpha*psi00 - y(3)*(y(1)+l)/y(2)+y(4))/(4*psi00);
u=0;
if (G > u_max^3)
    u = u_max;
else
    if (G < 0)
        u = 0;
    else
        u = G^(1/3);
    end
end
if (y(2)<=M)
    u = 0;
end

dydt(1) =-g+u*(y(1)+l)/y(2);
dydt(2) = -u;
dydt(3) = -y(3)*u/y(2)-psi30;
dydt(4) = y(3)*u*(y(1)+l)/(y(2)^2);
end

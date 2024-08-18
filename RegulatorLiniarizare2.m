%% Model Parameters
g = 9.80665;         % gravity constant
m_r = 0.2948;        % mass of the rod
m_w = 0.0695;        % mass of the inertia wheel
R = 0.05;            % radius of the inertia wheel
r = 0.02;            % cross section radius of the rod
l = 0.13;            % corresponding lengths
l_AD = l;
l_AC = l;            % assume wheel is mounted on the top of the pendulum
l_AB = l/2;
I_w_C = 0.5*m_w*R^2; % corresponding inertias
I_w_A = I_w_C + m_w*l_AC^2;
I_r_B = (1/12)*m_r*(3*r^2+l_AD^2);
I_r_A = I_r_B + m_r*l_AB^2;
Ts = 0.01;

c1 = (g*(m_r*l_AB + m_w*l_AC))/(I_w_A+I_r_A);
c2 = 1/(I_w_A+I_r_A);
c3 = 1/I_w_C;


%% State Space System -> Transfer Function

A = [0 1 0; 
     c1 0 0; 
     0 -1 0]
B = [0; -c2; c3]
C = [1 0 0] % we measure the angle of the pendulum rod - theta
D = 0

[num, den] = ss2tf(A,B,C,D);

H = minreal(zpk(tf(num, den)))


%%
t = 0:0.01:10;
u = t;
y = lsim(H,u,t);
plot(u,y)

%% Regulator cu G-T

Hf = tf(-338.8, [1, 0, -93.683])
sigma = 0.10; % suprareglaj 10%
tr = 3; % timp de raspuns 3 sec

zita = abs(log(sigma))/sqrt(log(sigma).^2 + pi^2)
wn = 4/(tr*zita)
cv = wn/(2*zita)
Estv = 1/cv
%wb = wn*sqrt(1 - 2*zita^2 + sqrt(2-4*zita^2+4*zita^4))

H0 = tf(wn^2, [1, 2*zita*wn, wn^2]);

t = 0: 0.1 : 10;
u = t;
lsim(u,t,H0)
figure
step(H0)
%figure;
%bode(H0);

Hr = 1/Hf * H0/(1-H0)
minreal(Hr)

%%

Te = 0.01;
Hfd = c2d(Hf, Te, 'zoh') % procesul discretizat
Hrd = c2d(minreal(Hr), Te, 'zoh') % regulatorul discretizat
H0d = feedback(series(Hrd,Hfd), 1) % H0 discretizat
step(H0d)
%%
t = 0:Te:10;
y = step(H0d,t)

plot(t,y), grid
%%
e = 1-y;
n = 267;
c(1) = -0.01502;
c(2) = -0.013703448;
c(3) = -0.0122482848;

%%
for k=4:n+1
    c(k) = 0.6703*c(k-1) - 0.01502*e(k) + 0.01502*e(k-1)
end

figure
stairs(t, c, 'r*')
hold on
plot(t,y),grid

%%
A = [0 1 0; 
     c1 0 0; 
     0 -1 0];
B = [0; -c2; c3];
C = [1 0 0]; % we measure the angle of the pendulum rod - theta
D = 0;

t = 0:0.01:10;
% Studiul controlabilitatii:
Co = ctrb(A,B);
rank(Co); % => rangul este 3, deci avem un sistem complet controlabil

% ne intereseaza Kp pentru pozitia pendulului si Kd pentru viteza acestuia 

sys_d = ss(A,B,C,D);
step(sys_d,t)


polii_impusi = [-1; -9.67; -9];
Kx = place(A, B, polii_impusi)

%Kx = [-10 -1 -0.0015];
sys_o = ss(A-B*Kx,B,C,D);
figure
step(sys_o,t),grid;

%% proiectare emipirica


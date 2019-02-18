% ONLY LATERAL PRESSURE MODEL NON LINEAR + HIGH PASS FILTER (percived pressure)+LOW PASS FILTER (valves dynamics)
% CONTACT RELATIONS
%( 6 STATI)
% ricordarsi di modificare anche InitData inizializzando opportunamente il
% vettore

%%***** SETTING MAIN MATMPC PATH

path_main_matmpc = 'C:\Users\auto\Documents\MATLAB\MATMPC';

%% Load params
%Pressure model params current

A = 0.016;          %area di contatto       
MM = 50;            %massa busto e braccia
m = 67;             %massa che interviene nella dinamica

k1 = 10*10^6; %12*10^5 %12*10^5;     %parametri forza elastica k1*(*dy)^p+k2  
k2 =100;
c1 = 4*10^6; %20*10^4;       %parametri smorzamento  c1*(*dy)^p+c2      
c2 = 800; %1000;
p=2;                %grado del polinomio smorzamento e rigidezza

alpha = 10;         %inclinazione poggiaschiena sedile
sigma_0 = 10^3;     %parametro per l'attrido [de Wit]
vs = 0.005;         %velocita` di Stribeck 
Fs = 40;            %coefficente di attrito statico
Fc = 30;            %coefficente di attrito dinamico (coulomb)
g = 9.81;



%% ****** FILTERS PARAMs *****

% HIGH PASS FILTER

tau_hp = 21; %[sec]
G_hp_ref = 1;       %hp gain to create the reference
G_hp_mpc = 1; %7;     %guadagno modello mpc percezione
% 
G_hp_mpc_a = 1; % 1/5.5; %guadagno modello mpc percezione active seat ( _v2)
G_hp_mpc_p = 1;% 7; %guadagno modello mpc percezione piattaforma ( _v2)

% LOW PASS FILTER

tau_lp = 0.1;
G_lp_mpc=1;


%%******* CONTACT PARAMs *******

% caratteristiche elastiche tronco(1) e cuscinetto(2)

E1 = 9.5*10^5;     %modulo di young [GPa] 
E2 = 0.1*10^5;
ni1 = 0.20;         %coefficente di poisson
ni2 = 0.48;

E = (E1/(1-ni1^2))+(E2/(1-ni2^2));      %E*=parametro che tiene in considerazione le caratteristiche elastiche dei corpi

% linearizzazione funzione R(u) u= pressione cuscinetto
q_lin = 10.8*10^5;
c_ang = 270*10^5;

%% Dimensions

nx=6;       % No. of states
nu=1;       % No. of controls
ny=2;       % No. of outputs
nyN=1;      % No. of outputs at the terminal point
np=3;       % No. of model parameters
nc=0;       % No. of general constraints
ncN=0;      % No. of general constraints at the terminal point
nbx = 0;    % No. of bounds on states
nbu = 0;    % No. of bounds on controls


import casadi.*

states   = SX.sym('states',nx,1);
controls = SX.sym('controls',nu,1);
params   = SX.sym('params',np,1);    
refs     = SX.sym('refs',ny,1);
refN     = SX.sym('refs',nyN,1);
Q        = SX.sym('Q',ny,ny);
QN       = SX.sym('QN',nyN,nyN);


%% Dynamics

accX=params(1); 
roll=params(2); 
accY=params(3);

prY1=states(1); 
prY2=states(2); 
prY3=states(3); 
pressY=states(4);
x_hp = states(5);
x_lp = states (6);

dpressY=controls(1);

   
%% versione 5 stati [no y_press_hp no y_press], espressione uscita su funzione di costo

tmp1= (sqrt(prY2^2)*prY3) ; %modulo( prY2 )* prY3
tmp2= m*accX*cos(pi/180*alpha)+MM*g*sin(pi/180*alpha) ; %Fn = normale al poggia schiena
tmp3= 1/(pi)*atan(tmp2)+0.6 ;   %fuzione per annullare o meno attrito in presenza di contatto o meno con il sedile

tmp4 = -(c1*(prY1)^2+c2)/m*prY2-(k1*(prY1)^2+k2)*prY1/m+accY+MM*g*roll/m-sigma_0*prY3/m; %body dynamic eq

x_dot=[prY2;...
       tmp4; ...
       prY2-tmp1/((Fc*tmp3+((Fs-Fc)*tmp3*exp(-(prY2/vs)^2)))/sigma_0);...               
       dpressY;... %per non avere il controllo in uscita
       (-1/tau_hp)*x_hp+[(c1*(prY1)^2)*prY2+(k1*(prY1)^2)*prY1]/A+[(1/tau_lp)*G_lp_mpc*x_lp];...
       (-1/tau_lp)*x_lp+pressY];

% x_dot=[prY2;...
%        tmp4; ...
%        prY2-tmp1/((Fc*tmp3+((Fs-Fc)*tmp3*exp(-(prY2/vs)^2)))/sigma_0);...               
%        (-1/tau_lp)*x_lp+dpressY;... %per non avere il controllo in uscita
%        (-1/tau_hp)*x_hp+[(c1*(prY1)^2)*prY2+(k1*(prY1)^2)*prY1]/A+[(1/tau_lp)*x_lp]];...
   
%%

xdot = SX.sym('xdot',nx,1);
impl_f = xdot - x_dot;
   
 
xdot = SX.sym('xdot',nx,1);
impl_f = xdot - x_dot;
     
%% Objectives and constraints

% objectives

h = [(-1/tau_hp)*G_hp_mpc*x_hp+G_hp_mpc*([(c1*(prY1)^2)*prY2+(k1*(prY1)^2)*prY1]/A+[(1/tau_lp)*G_lp_mpc*x_lp]); pressY];

hN=[pressY]; %generic state 

h_fun=Function('h_fun', {states,controls,params}, {h},{'states','controls','params'},{'h'});
hN_fun=Function('hN_fun', {states,params}, {hN},{'states','params'},{'hN'});

% general inequality constraints
general_con = [];
general_con_N = [];

% state and control bounds
nbx_idx = 0;    %nbx_idx = [2]; % indexs of states which are bounded
nbu_idx = 0;    % indexs of controls which are bounded
path_con=general_con;
path_con_N=general_con_N;
for i=1:nbx
    path_con=[path_con;states(nbx_idx(i))];
    path_con_N=[path_con_N;states(nbx_idx(i))];
end    
nc=nc+nbx;
ncN=ncN+nbx;

% build the function for inequality constraints
path_con_fun=Function('path_con_fun', {states,controls,params}, {path_con},{'states','controls','params'},{'path_con'});
path_con_N_fun=Function('path_con_N_fun', {states,params}, {path_con_N},{'states','params'},{'path_con_N'});


%% NMPC sampling time [s]

Ts = 0.005; % simulation sample time
Ts_st = 0.005; % shooting interval time

%% save your data in the path of your MATMPC

cd(path_main_matmpc); %main matmpc folder

clc;

% %% Dimensions
% 
% nx=6;       % No. of states
% nu=1;       % No. of controls
% ny=2;       % No. of outputs
% nyN=1;      % No. of outputs at the terminal point
% np=3;       % No. of model parameters
% nc=0;       % No. of general constraints
% ncN=0;      % No. of general constraints at the terminal point
% nbx = 0;    % No. of bounds on states
% nbu = 0;    % No. of bounds on controls
% 
% 
% import casadi.*
% 
% states   = SX.sym('states',nx,1);
% controls = SX.sym('controls',nu,1);
% params   = SX.sym('params',np,1);    
% refs     = SX.sym('refs',ny,1);
% refN     = SX.sym('refs',nyN,1);
% Q        = SX.sym('Q',ny,ny);
% QN       = SX.sym('QN',nyN,nyN);
% 
% 
% %% Dynamics
% 
% accX=params(1); 
% roll=params(2); 
% accY=params(3);
% 
% prY1=states(1); 
% prY2=states(2); 
% prY3=states(3); 
% pressY=states(4);
% x_hp=states(5);
% x_lp=states(6);
% 
% dpressY=controls(1);
% 
%    
% %% versione 6 stati 
% 
% tmp1= (sqrt(prY2^2)*prY3) ; %modulo( prY2 )* prY3
% tmp2= m*accX*cos(pi/180*alpha)+MM*g*sin(pi/180*alpha) ; %Fn = normale al poggia schiena
% tmp3= 1/(pi)*atan(tmp2)+0.6 ;   %fuzione per annullare o meno attrito in presenza di contatto o meno con il sedile
% 
% tmp4 = -(c1*(prY1)^2+c2)/m*prY2-(k1*(prY1)^2+k2)*prY1/m+accY+MM*g*roll/m-sigma_0*prY3/m; %body dynamic eq
% 
% 
% % temp5 = dpressY/c_ang;                      %d=profondita` espressa come la derivata del raggio di curvatura
% % temp6 = (pressY+q_lin)/c_ang;               %R= raggio di curvatura in funzione della pressione
% % temp7 = (2/pi)*E*(temp5/temp6)^(1/2);       %Ps= pressione subita dal corpo nel contatto (formula di Hertz)
% 
% x_dot=[prY2;...
%        tmp4; ...
%        prY2-tmp1/((Fc*tmp3+((Fs-Fc)*tmp3*exp(-(prY2/vs)^2)))/sigma_0);...               
%        dpressY;... %per non avere il controllo in uscita
%        (-1/tau_hp)*x_hp+[(c1*(prY1)^2)*prY2+(k1*(prY1)^2)*prY1]/A+[(1/tau_lp)*G_lp_mpc*x_lp];...
%        (-1/tau_lp)*x_lp+pressY];
% 
% % x_dot=[prY2;...
% %        tmp4; ...
% %        prY2-tmp1/((Fc*tmp3+((Fs-Fc)*tmp3*exp(-(prY2/vs)^2)))/sigma_0);...               
% %        dpressY;... 
% %        (-1/tau_hp)*x_hp+((c1*(prY1)^2)*prY2+(k1*(prY1)^2)*prY1)/A+((1/tau_lp)*G_lp_mpc*x_lp);...
% %        (-1/tau_lp)*x_lp+(2/pi)*E*((prY1-(dpressY/c_ang))/((pressY+q_lin)/c_ang))^(1/2)];
% 
%    
% %%
%  
% xdot = SX.sym('xdot',nx,1);
% impl_f = xdot - x_dot;
%      
% %% Objectives and constraints
% 
% % objectives
% h = [(-1/tau_hp)*G_hp_mpc*x_hp+G_hp_mpc*([(c1*(prY1)^2)*prY2+(k1*(prY1)^2)*prY1]/A+[(1/tau_lp)*G_lp_mpc*x_lp]); pressY];
% 
% %h = [(-1/tau_hp)*x_hp+[(c1*(prY1)^2)*prY2+(k1*(prY1)^2)*prY1]/A+ pressY; pressY];
% %(-1/tau_hp)*G_hp_mpc*x_hp+G_hp_mpc*(((c1*(prY1)^2)*prY2+(k1*(prY1)^2)*prY1)/A+((1/tau_lp)*G_lp_mpc*x_lp))
% hN=[pressY]; %generic state 
% 
% h_fun=Function('h_fun', {states,controls,params}, {h},{'states','controls','params'},{'h'});
% hN_fun=Function('hN_fun', {states,params}, {hN},{'states','params'},{'hN'});
% 
% % general inequality constraints
% general_con = [];
% general_con_N = [];
% 
% % state and control bounds
% nbx_idx = 0;    %nbx_idx = [2]; % indexs of states which are bounded
% nbu_idx = 0;    % indexs of controls which are bounded
% path_con=general_con;
% path_con_N=general_con_N;
% for i=1:nbx
%     path_con=[path_con;states(nbx_idx(i))];
%     path_con_N=[path_con_N;states(nbx_idx(i))];
% end    
% nc=nc+nbx;
% ncN=ncN+nbx;
% 
% % build the function for inequality constraints
% path_con_fun=Function('path_con_fun', {states,controls,params}, {path_con},{'states','controls','params'},{'path_con'});
% path_con_N_fun=Function('path_con_N_fun', {states,params}, {path_con_N},{'states','params'},{'path_con_N'});
% 
% 
% %% NMPC sampling time [s]
% 
% Ts = 0.005; % simulation sample time
% Ts_st = 0.005; % shooting interval time
% 
% %% save your data in the path of your MATMPC
% 
% cd(path_main_matmpc); %main matmpc folder
% 
% clc;
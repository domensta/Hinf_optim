

close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcule un correcteur Hinfini � partir du
% sch�ma bloc d�fini dans master_2005.mdl
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% chargement des pond�rations
pond = [];
save pond;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% r�cup�ration du mod�le simulink
[a,b,c,d]=linmod('modele_syn');
close all

% mod�le de synth�se
[a1,b1,c1,d1]=minreal(a,b(:,1:3),c(1:3,:),d(1:3,1:3));
ssaug=ss(a1,b1,c1,d1);
sysaug=ltisys(a1,b1,c1,d1);

load('pond.mat', 'inv_W1');
W1 = 1/inv_W1;
load('pond.mat', 'inv_W2');
W2 = 1/inv_W2;
load('pond.mat', 'inv_W3');
W3 = 1/inv_W3;
s = tf('s');
Ir = 40;
ws = 0.06
zs = 0.0005;
L = -2.5;
G = tf([1 2*ws*zs ws*ws],[Ir-L*L Ir*2*ws*zs Ir*ws*ws 0]);
G = G/s;


%% hinfstruct synthesis


% define PID controller, 'P', 'PI', 'PD' also possible
 options = hinfstructOptions('RandomStart',10, 'Display', 'off');
C0pid = ltiblock.pid('namePID', 'PID') ; %
%C0pid.Tf.Free=0;
%C0pid.Tf.Value=1;
C0pid.Tf.Minimum=0;
C0pid.Tf.Maximum=5;
C0pid.Kp.Minimum=-10;
C0pid.Ki.Minimum=-10;
C0pid.Kd.Minimum=-10;
C0pid.Kp.Maximum=10
C0pid.Ki.Maximum=10;
C0pid.Kd.Maximum=10;

% Cpid is tuned version of C0pid, gam is optimal Hinf norm
[Cpid, GAM, info] = hinfstruct(ssaug, C0pid,options) ;
trace_hinf(ssaug,Cpid,inv_W2,inv_W3,inv_W1);
p = pole(Cpid);
if sum(real(p)>0)>0,
   disp('CORRECTEUR INSTABLE !!!!!')
end;
%% hinf synthesis
[syscor,gopt]=hinfsyn(sysaug,1,1,0.9,10,0.1);

[acor,bcor,ccor,dcor]=ltiss(syscor);
cor=ss(acor,bcor,ccor,dcor);
[num_cor,den_cor]=tfdata(cor,'v');

p=pole(cor);
if sum(real(p)>0)>0,
   disp('CORRECTEUR INSTABLE !!!!!')
end;

% v�rification des contraintes
% sys : 	entr�es = [ consigne ; Cpert ; cor]
%			sorties = [teta ; erreur ]
% cor : 	entr�e =  [erreur]
%	  :	    sortie = sys(3)
[a2,b2,c2,d2]=minreal(a,b(:,[1 4 3]),c([4 3],:),d([4 3],[1 4 3]));
sys=ss(a2,b2,c2,d2);

% perfo Hinfini

trace_hinf(sys,cor,inv_W2,inv_W3,inv_W1)

pid =4.93953 + 0.0490246/s +7.16762*s/(1+s); % not remove
0.975822 ; 0.000988356 ; 9.41187 ; 3.78795
pid =1.0485 + 0.0018779/s +9.49405*s/(1+3.29043*s);
trace_hinf(sys,pid,inv_W2,inv_W3,inv_W1);
trace_hinf(sys,Cpid,inv_W2,inv_W3,inv_W1);
% marges de stabilit�
figure
nichols(sys(1,3)*pid)
hold on
plot(-180,0,'k+')
grid
zoom on
hold off
%%
mat = [-88.8577 , -1973.42 , 0.223691 , 0 , 0 , 0 , 0 , 0 , 4.70394 , 0.00142706;
-1973.42 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , -0.5;
0.223691 , 0 , 0 , -0.5 , 0 , 0.238506 , 8.88889e-07 , 5.33333e-05 , 0 , 0;
0 , 0 , -0.5 , -0.0280993 , 0.499803 , 0 , 0 , 0 , 0 , 0;
0 , 0 , 0 , 0.499803 , 0 , 0 , 0.5 , 0 , 0 , 0;
0 , 0 , 8.88889e-07 , 0 , 0.5 , 0 , 0 , 0 , 0 , 0;
0 , 0 , 5.33333e-05 , 0 , 0 , 0.5 , 0 , 0 , 0 , 0;
4.70394 , 0 , 0 , 0 , 0 , 4.70394 , 0 , 0 , 1 , 0;
0 , 0 , 0.238506 , 0 , 0 , 0 , 0 , 0.5 , 4.70394 , 0.00142706;
0.00142706 , -0.5 , 0 , 0 , 0 , 0.00142706 , 0 , 0 , 0 , -0.994219]

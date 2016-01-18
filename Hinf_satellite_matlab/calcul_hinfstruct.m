

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


%% hinfstruct synthesis

% define PID controller, 'P', 'PI', 'PD' also possible
 options = hinfstructOptions('RandomStart',1, 'Display', 'off');
C0pid = ltiblock.pid('namePID', 'PID') ; %
C0pid.Tf.Free=0;
C0pid.Tf.Value=1;
C0pid.Kp.Minimum=-10;
C0pid.Ki.Minimum=-10;
C0pid.Kd.Minimum=-10;
C0pid.Kp.Maximum=10;
C0pid.Ki.Maximum=10;
C0pid.Kd.Maximum=10;

% Cpid is tuned version of C0pid, gam is optimal Hinf norm
[Cpid, GAM, info] = hinfstruct(ssaug, C0pid,options) ;
norm(lft(Cpid,ssaug,1,1),inf);

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

% marges de stabilit�
figure
nichols(sys(1,3)*cor)
hold on
plot(-180,0,'k+')
grid
zoom on
hold off

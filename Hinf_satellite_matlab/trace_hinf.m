function trace_hinf(sys,cor,inv_W2,inv_W3,inv_W1)

u = -3:0.0001:3;
w = 10.^u;
% sys : 	entr�es = [ consigne ; Cpert ; cor]
%			sorties = [teta ; erreur  ]
% cor : 	entr�e = [erreur]
%		:	sortie = sys(3)

% construction du mod�le boucl� :
% entr�es : [consigne, Cpert]
% sorties : [teta ; erreur ; cor]


% concat�nation tous syst�me boucle ouverte
syscat=append(sys,cor);

% matrice d'interconnection
%	entr�e (1) sys = pendante (0)
%	entr�e (2) sys = pendante (0)
%	entr�e (3) sys = cor (3)
%	entr�e (4) cor = erreur (2)

Q=[1 0;2 0;3 3;4 2];


% entr�es globales du syst�me [consigne;cpert]
in=[1 2];


% sorties globales du syst�me [teta ; erreur ; cor]
out=[1 2 3];

% syst�me augment�
sysbf=connect(syscat,Q,in,out);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRACES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure

% erreur/teta_cons compar� � inv_W1
subplot(2,2,1)
[mag,ph]=bode(sysbf(2,1),w);
[mage,phe]=bode(inv_W1,w);
semilogx(squeeze(w),20*log10(squeeze(mag)),squeeze(w),20*log10(squeeze(mage)))
title('e / \theta_c, inv(W_1)')
zoom on
grid

% cor/cons compar� � inv_W2
subplot(2,2,2)
[mag,ph]=bode(sysbf(3,1),w);
[magk,phk]=bode(inv_W2,w);
semilogx(squeeze(w),20*log10(squeeze(mag)),squeeze(w),20*log10(squeeze(magk)))
title('u / \theta_c, inv(W_2)')
zoom on
grid

% teta/cpert compar� � inv_W1*inv_W3
subplot(2,2,3)
[mag,ph]=bode(sysbf(1,2),w);
[mageg,pheg]=bode(inv_W1*inv_W3,w);
semilogx(squeeze(w),20*log10(squeeze(mag)),squeeze(w),20*log10(squeeze(mageg)))
title('e / Cpert, inv(W_1*W_3)')
zoom on
grid

% teta/consigne compar� � inv_W2*inv_W3
subplot(2,2,4)
[mag,ph]=bode(sysbf(1,1),w);
[magkg,phkg]=bode(inv_W2*inv_W3,w);
semilogx(squeeze(w),20*log10(squeeze(mag)),squeeze(w),20*log10(squeeze(magkg)))
title('u/Cpert = \theta / \theta_c, inv(W_2*W_3)')
zoom on
grid



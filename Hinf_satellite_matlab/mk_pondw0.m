function px=mk_pondw0(gainbf,coupure,pente,gainhf,titre)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% px=mk_pond(gainbf,coupure,pente,gainhf)
%
% fonction fabriquant les pondérations H infini monovariables 
% à partir de filtres élémentaires
%
% gainbf = vecteur des gains en basse fréq (en dB)
% coupure = vecteur des fréq de passage à 0 dB(Hz)
% pente = vecteur des pentes (= ordre des filtres)
% gainhf = vecteur des gains en haute fréq (en dB)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% nombre de filtres  
n_filtre = length(gainbf);
px=ss(1);

if nargin<5,
    titre='';
end;

for ii=1:n_filtre,
   
      
   % filtre passe-bas
   if pente(ii) <0,
      
      ordre=abs(pente(ii));
      d_gain=gainhf(ii)-gainbf(ii);
      coupure_bf(ii)=10^(gainbf(ii)/20/pente(ii))*coupure(ii);
      [n_pole,d_pole]=butter(ordre,coupure_bf(ii)*2*pi,'s');
      [n_zero,d_zero]=butter(ordre,10^(d_gain/20/pente(ii))*coupure_bf(ii)*2*pi,'s');
      

	  if ii>1,
   		 num=d_zero/d_zero(ordre+1);
	  else
         num=10^(gainbf(ii)/20)*d_zero/d_zero(ordre+1);
      end;
      den=d_pole/d_pole(ordre+1);
      
   else
   % filtre passe-haut
   
   	ordre=abs(pente(ii));
      d_gain=gainhf(ii)-gainbf(ii);
      coupure_bf(ii)=10^(gainbf(ii)/20/pente(ii))*coupure(ii);
      [n_pole,d_pole]=butter(ordre,coupure_bf(ii)*2*pi,'s');
      [n_zero,d_zero]=butter(ordre,10^(d_gain/20/pente(ii))*coupure_bf(ii)*2*pi,'s');
      
      if ii>1,
         num=d_pole/d_pole(ordre+1);
      else
         num=10^(gainbf(ii)/20)*d_pole/d_pole(ordre+1);
		end;
      den=d_zero/d_zero(ordre+1);
      
      
	end;
   
   pxe=tf(num,den); 
   px=series(px,pxe);

end;

[mag,ph,w]=bode(px);
subplot(2,1,1)
semilogx(w/2/pi,20*log10(squeeze(mag)));
ylabel('dB')
title(['Pondération ',titre])
grid
subplot(2,1,2)
semilogx(w/2/pi,squeeze(ph));
grid
xlabel('Hz')
ylabel('deg')
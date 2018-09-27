function [der_S,der_I,der_E] = VecFieldNoParODE (S, I, E, tau)
  %Funzione per il calcolo del campo vettoriale associato al modello SIR con vaccinazione
  
  n=length(E);
  
  %--------Parametri e funzioni--------%
  
  mu=1/(75*365);         %giorni^(-1), l'inverso della speranza di vita
  p0=0.75;               %baseline
  nu=1/7;                %giorni^(-1), l'inverso della durata media della malattia
  beta=20*(mu+nu);       %R0=20;
  a=n/tau;               %parametro della distribuzione di Erlang
  c=300;                 %parametro di p1;
  p1=min(c*E(n),1-p0);   
  k=1;                   %parametro di g;
  g=k*I;                 
  
  %--------Calcolo del campo vettoriale--------%
  
  der_S=mu*(1-p0-p1)-mu*S-beta*S*I;
  der_I=I*(beta*S-mu-nu);
  der_E(1)=a*g-a*E(1);
  for i=2:n
    der_E(i)=a*E(i-1)-a*E(i);
  end
  
end

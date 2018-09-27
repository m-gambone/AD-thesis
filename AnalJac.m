function [Jac] = AnalJac (S, I, E, tau)
  %Funzione per il calcolo analitico della Jacobiana del modello SIR con vaccinazione
  
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
  
  Jac(1,:)=[-mu-beta*I, -beta*S, zeros(1,n-1), -mu*c];
  Jac(2,:)=[beta*I, beta*S-(mu+nu), zeros(1,n)];
  Jac(3,1:2)=[0,k*a];
  Jac(3:n+2,3:n+2)=diag(-a*ones(n,1),0)+diag(a*ones(n-1,1),-1);
  
 end
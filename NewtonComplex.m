function [alfa,k]=NewtonComplex(x0)
  %Metodo di Newton per calcolare gli equilibri del modello SIR con vaccinazione (admDiffComplex)
  
  %--------Parametri e variabili--------%
  
  TOL=eps;
  kmax=100;
  x=x0;
  k=1;
  error(1)=TOL+1;
  y=x0;
  
  %--------Opzioni di ADiMat--------%
  
  adopts = admOptions('i',[1:3]);
  adopts.flags = '--check-certificate';
  
  %--------Calcolo dell'equilibrio--------%
  
  while (k<kmax) & (error>=TOL)
    S=x(1);
    I=x(2);
    E=x(3:length(x));
    Jac=admDiffComplex(@VecFieldNoParODE, 1, S, I, E, 1, adopts);
    [y(1),y(2),y(3:length(y))]=VecFieldNoParODE(S, I, E, 1);
    dx=(-Jac\y')';
    x=x+dx;
    error(k+1)=norm(dx);
    k=k+1;
  end
  alfa=x;
  
end  
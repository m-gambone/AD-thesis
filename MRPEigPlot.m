function [y]=MRPEigPlot()
  %Funzione per il calcolo del max(real(eig(Jac))) e produzione relativo grafico
  
  %--------Parametri e variabili--------%
  
  ParDiscr=50;
  x=linspace(30,10000,ParDiscr);
  
  %--------Opzioni di ADiMat--------%
  
  adopts = admOptions('i',[1:3]);
  adopts.flags = '--check-certificate';
  
  %--------Calcolo del max(real(eig(Jac)))--------%
  
  for n=[1:5]
    equil=NewtonVFor(-1*ones(1,n+2));
    for tau=1:ParDiscr
      S=equil(1);
      I=equil(2);
      E=equil(3:length(equil));
      Jac=admDiffVFor(@VecFieldNoParODE, 1, S, I, E, x(tau), adopts);
      y(n,tau)=max(real(eig(Jac)));
    end
    %--------Grafico max(real(eig(Jac)))--------%  
    plot(x,y(n,:),'-')
    hold on
  end
  
  plot(x,zeros(1,length(x)),'--')
  xlabel('$$\tau$$','Interpreter','latex')
  ylabel('$$\max_i\Re(\lambda_i)$$','Interpreter','latex')
  legend({'n=200','n=400','n=800'},'Location','southeast')
  
end
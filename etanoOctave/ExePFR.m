clear all; close all; clc;
warning off;
global A11 A12 A2 A3 A41 A42 E11 E12 E2 E3 E41 E42
global RG1 RG2
global T P Qf

#Prâmetros Cinéticos
A11=1e14; A12=1e12; A2=3e14; A3=3.4e12; A41=1e12; A42=1e13;
E11=217.6e3; E12=0; E2=165.3e3; E3=28.5e3; E41=0; E42=200.8e3;
#Constantes Físicas
RG1=8.314; # j/mol K
RG2=82.06; #cm³ atm/ mol K
#Constantes de operação
T = 1050;  #K
P = 1;     #atm 
Qf = 600;  #cm^3/s

#funÃ§Ã£o para definiÃ§Ã£o das equaÃ§Ãµes diferenciais a serem integradas
function dxdt = fun(y,t)
  global A11 A12 A2 A3 A41 A42 E11 E12 E2 E3 E41 E42
  global RG1 RG2
  global T P Qf
  #t ==> variÃ¡vel independent
  #y ==> vetor de variÃ¡veis dependentes
  #dxdt ==> vetor dos valores das EDOs em t e y
  #Área de identificação de variáveis
      #legenda
      # N1 --> C2H6
      # N2 --> C2H5⁺
      # N3 --> C2H4
      # N4 --> H⁺
      # N5 --> H2
      # N6 --> NO
      # N7 --> HNO 
  V=t;
  N1=y(1); N2=y(2); N3=y(3); N4=y(4);
  N5=y(5); N6=y(6); N7=y(7);
  #Área para definição de parâmtros

  #Equações auxiliares
  k11=A11*exp(-E11/(RG1*T)); k12=A12*exp(-E12/(RG1*T));
  k2=A2*exp(-E2/(RG1*T)); 
  k3=A3*exp(-E3/(RG1*T));
  k41=A41*exp(-E41/(RG1*T)); k42=A42*exp(-E42/(RG1*T));  
  
  Q=(RG2*T/P)*(N1+N2+N3+N4+N5+N6+N7);
  
  C1=N1/Q; C2=N2/Q; C3=N3/Q; C4=N4/Q;
  C5=N5/Q; C6=N6/Q; C7=N7/Q;
  
  r1=k11*C1*C6-k12*C2*C7;
  r2=k2*C2;
  r3=k3*C4*C1;
  r4=k41*C4*C6-k42*C7;
  
  R1=-r1-r3; R2=r1-r2+r3; R3=r2; R4=r2-r3-r4;
  R5=r3; R6=-r1-r4; R7=r1+r4;
  
  #DefiniÃ§Ã£o das EDOs
  dxdt(1)=R1;
  dxdt(2)=R2;
  dxdt(3)=R3;
  dxdt(4)=R4;
  dxdt(5)=R5;
  dxdt(6)=R6;
  dxdt(7)=R7;
endfunction


#definiÃ§Ã£o de parÃ¢metros, condiÃ§Ãµes iniciais e intervalo de integraÃ§Ã£o
#intervalo de integraÃ§Ã£o
t0=0;
tf=1500;
trange=t=[t0:tf/100:tf];;
#CondiÃ§Ãµes inicias
y0=[0.95*Qf*P/(RG2*T) 0 0 0 0 0.05*Qf*P/(RG2*T) 0];
#parÃ¢metros numÃ©ricos
y = lsode(@fun,y0,trange);

Q = (RG2*T/P)*(sum(y,2));

#AnÃ¡lise de dados
#identificaÃ§Ã£o de variÃ¡veis
V=trange;
N1=y(:,1);
N2=y(:,2);
N3=y(:,3);
N4=y(:,4);
N5=y(:,5);
N6=y(:,6);
N7=y(:,7);

#Normalização
N_Normalization = y./max(y);

figure(1)
plot(V,N1,'-b;;','linewidth',3,...
     V,N3,'-r;;','linewidth',3,...
     V,N6,'-g;;','linewidth',3)
h=legend({"C_2H_6","C_2H_4","N_O"});
set(h,'fontsize',14)
legend boxoff;
set(gca,'fontsize',12)
grid on;
xlabel("Volume [cm³]");
ylabel("N_j [mol/s]")
axis([t(1) t(end)]);


figure(2)
plot(V,N_Normalization,'linewidth',3)
h=legend({"C_2H_6","C_2H_5^+","C_2H_4","H^+","H_2","N_O","HNO"});
legend boxoff;
set(gca,'fontsize',12)
set(h,'fontsize',14,'location','northeastoutside');
grid on;
xlabel("Volume [cm³]");
ylabel("N_j/N_{j max} [mol/s]")
axis([t(1) t(end)]);

figure(3)
plot(V,Q,'linewidth',3)
set(gca,'fontsize',12)
grid on;
xlabel("Volume [cm³]");
ylabel("Taxa Volumetrica[cm³/s]")
axis([t(1) t(end)]);


#Variação de Temperatura
t1050 = trange;
N11050 = N1;

clear y;
T=1100;
y=lsode(@fun,y0,trange);
t1100 = trange;
N11100 = y(:,1);

clear y;
T=1000;
y=lsode(@fun,y0,trange);
t1000 = trange;
N11000 = y(:,1);

figure(4)
plot(t1000,N11000,'-b;;','linewidth',3,...
     t1050,N11050,'-r;;','linewidth',3,...
     t1100,N11100,'-g;;','linewidth',3)
h=legend({"T = 1000 K","T = 1050 K","T = 1100 K"});
set(h,'fontsize',14)
legend boxoff;
set(gca,'fontsize',12)
grid on;
xlabel("Volume [cm³]");
ylabel("N_{C_2H_6} [mol/s]")
axis([t(1) t(end)]);



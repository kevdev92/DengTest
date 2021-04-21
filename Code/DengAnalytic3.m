%%
clear all
close all

%%%%% Comparing data from figures in Deng et al, analytic solution given by
%%%%% Deng et al and COMSOL solution to stress PDE given in Deng et al.

format long
set(0,'defaulttextinterpreter','Latex')

tic
%%%Paramaters
RO = 4e-8;
delta = 0.03; %SEI width
M = 1001; %space
x = 1e-16/RO:2e-10/RO:RO/RO; %active region 
x2 = 1:8e-4:1+delta;
N = 4; %time
k = 100; % number of terms in bessel series, the larger the number, the longer the code will take to run

%%%% Import data from Deng paper, 
M6 = readtable('Deng3aData.csv'); %Deng Data stress
xx6 = table2array(M6(1:end,1)); %x data
sigDeng = table2array(M6(1:end,2)); %stress data

M8 = readtable('Deng2Data.csv'); %Deng Data Concentration
xx8 = table2array(M8(1:end,1)); % x data
CCDengDat = table2array(M8(1:end,2)); %concentration data

%%% Import Data from COMSOL computations
M4 = readtable('SEILin_PR_Om_Zero.csv'); %Comsol Linear model data
xx4 = table2array(M4(8:end,1)); %x data
uu4 = table2array(M4(8:end,2)); %stress data

nua = 0.2; %poissons ratio active material
nus = 0.3; %poissons ratio SEI
Ea = 1e10; %Young's modulous active material
ESEI = 1e9; %Young's modulous active material
beta = (ESEI/Ea)^1; %ratio of Young's modulii 
Omegaa = 8.9e-6; %partial molar volume active material
OmegaSEI = 0*8.9e-6; %partial molar volume SEI
alpha = OmegaSEI/Omegaa; %ratio between partial molar volumes

Sigmasur = -0.001/3/(1-nua); %surface stress
Sigmasur = 0; %surface stress
D = 1e-14; %diffusion coefficient

Jc = 0.01; %current flux density
t = linspace(0,0.2,5); %time vector
t = [0.05 0.1 0.2];
C = zeros(length(x),length(t)); %initialise Concentration in active material
SSEI = zeros(length(x2),length(t));
Sa = C; %initialise Stress active material
CSEI = zeros(1,length(t)); %initialise Concentration in active material

%%% calculate analytic result using bessel functions
for j = 1:length(t)
    for i = 1:length(x) %%active material
        C(i,j) = (0.5*x(i)^2+2*(t(j)) - 0.25 - 2*sumbessel3(x(i),k,t(j))); 
        Sa(i,j) = Sigmaafun(nus, nua, Sigmasur,beta,alpha,delta,k,t(j),x(i)); %function to determine stress in active material
    end
    CSEI(j) = C(end,j);
    for i = 1:length(x2) %SEI
        SSEI(i,j) = SigmaSEIfun(nus, nua, Sigmasur,beta,alpha,delta,k,t(j),x2(i)); %function to determine stress in SEI

    end
end
toc
%%
tic
close all

%%%%%concnetration plots
ii = 1; % initialise counter to be used in loop for the legend
leg = 0; % initialise legend string variable
figure(1)
hold on
for j = 1:length(t)
plot([x x2],[C(:,j) ; CSEI(j)*ones(length(x2),1)],'Linewidth',1.2)
leg(ii) = t(j); %set legend string variable
ii = ii+1;
end
plot(xx8,CCDengDat,'kx')
ylim([-0.1 0.7])
xlim([0 1.03])
xlabel('$r/R$ ')
ylabel('$C/(RJ_C/D)$')
legendStrings = ["$Dt/R^{2} = $ " + string(leg) + " s", "Data ($Dt/R^{2} = 0.2$ s)"]; %legend adapts to number of time steps
legend(legendStrings,'Interpreter','latex','FontSize',11,'Location','northwest')

%%% plot comparing Deng plot, with COMSOL solution and analytic solution
%%% given in Deng
figure(2)
hold on
plot(xx6,sigDeng,'xk','Linewidth',1)
plot([x x2],[Sa(:,end);SSEI(:,end)]*3*(1-nua),'Linewidth',1.2,'color',[0 0.447 0.741])
plot(xx4/RO,uu4,'-.','Linewidth',1.2,'color',[0.85 0.325 0.098])
%grid on
%ylim([-0.1 0.7])
xlim([0 1+delta])
xlabel('$r/R$ ')
ylabel('$3(1-\nu_a)\Sigma_R$')
bleg2 = legend('Deng \textit{et al.} data', 'Deng \textit{et al.} analytic','Deng \textit{et al.} \texttt{COMSOL}'); 
set(bleg2,'Interpreter','Latex','Fontsize',11,'Location','northeast');
axes('Position',[.2 .2 .35 .35])
box on
hold on
plot(xx6,sigDeng,'xk','Linewidth',1)
plot([x x2],[Sa(:,end);SSEI(:,end)]*3*(1-nua),'Linewidth',1.2,'color',[0 0.447 0.741])
plot(xx4/RO,uu4,'-.','Linewidth',1.2,'color',[0.85 0.325 0.098])
grid on
ylim([-0.005 0.005])
xlim([0.975 1+delta])


toc
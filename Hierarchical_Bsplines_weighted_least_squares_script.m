clear all
close all
format long

%Definiamo il dominio su cui sono definite le B-spline bivariate: [0,1]x[0,1]
problem_data.geo_name = 'geo_square.txt';

% Qui inizia una parte di definizione di parametrici tecnici, necessari
% per usare geopdes per risolvere equazioni alle derivate parziali,
% ma che per noi non sono importanti, salvo grado e ordine di continuità
% e poche altre cose (si veda sotto)
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [1 2 3 4];

problem_data.c_diff  = @(x, y) ones(size(x));
C = 100;
normax2 = @(x,y) ((x-.5).^2+(y-.5).^2);
problem_data.uex = @(x,y) exp(-C*normax2(x,y));
problem_data.f = @(x,y) 4*C*(1-C*normax2(x,y)).*problem_data.uex(x,y);
problem_data.g = @(x, y, ind) zeros(size(x));
problem_data.h = @(x, y, ind) problem_data.uex(x,y);

problem_data.uex =@(x,y) exp(-C*normax2(x,y));
problem_data.graduex = @(x,y) -2*C*cat (1, ...
            reshape (problem_data.uex(x,y).*(x-.5), [1, size(x)]), ...
            reshape (problem_data.uex(x,y).*(y-.5), [1, size(x)]));
        
clear method_data
%Tutti questi parametri sono coppie, 
%perché sono definiti in direzione x e y
method_data.degree      = [2 2];   % grado delle B-spline
method_data.regularity  = [1 1];   % ordine di continuità
method_data.nsub_coarse = [32 32];  % numero di sotto-intervalli iniziali della partizione
method_data.nsub_refine = [2 2];  % in quante parti dividiamo ogni sotto-intervallo quando raffiniamo
method_data.nquad       = [3 3];  % non importante per i nostri scopi
method_data.space_type  = 'standard'; % usiamo sempre 'standard' 
method_data.truncated   = 1;   % 0: spline gerarchiche non troncate, 1: spline gerarchiche troncate

% ADAPTIVITY PARAMETERS
clear adaptivity_data
%adaptivity_data.flag = 'elements';
adaptivity_data.flag = 'functions';  %raffiniamo indicando i supporti (di funzioni) da raffinare
adaptivity_data.mark_param = .5;   % non importante per i nostri scopi
adaptivity_data.mark_strategy = 'MS'; % non importante per i nostri scopi
%Parametri che determinano l'arresto dell'algoritmo
adaptivity_data.max_level = 5;  %numero massimo di livelli nello spazio gerarchico
adaptivity_data.max_ndof = 50000;  %dimensione massimo dello spazio gerarchico
adaptivity_data.num_max_iter = 11;  %numero massimo di iterazioni
adaptivity_data.max_nel = 50000;  %numero massimo di elementi nella mesh gerarchica
adaptivity_data.tol = 30;  % tolleranza sull'errore
adaptivity_data.adm = 0;  % non importante per i nostri scopi (per ora)

% parametri per ottenere grafici alla fine dell'esecuzione dell'algoritmo
plot_hmesh = true;
plot_discrete_sol = false;

%Input di alcuni esempi di insiemi di dati
exm = input('Esempio # (1 = black forest, 2 = glacier):');
%I dati sono nella forma (u,v,f(u,v))
switch exm
    case 1 
        fw = fopen('black_forest.txt','r'); 
    case 2
        fw = fopen('glacier.txt','r');
end

data_3=fscanf(fw,'%f %f \r\n',[3 inf]);
fclose(fw);

data=data_3(1:2,:)';
%f=data_3(3,:)'/max(abs(data_3(3,:)));
f=data_3(3,:)';
a=min(data(:,1));
b=max(data(:,1));
c=min(data(:,2));
d=max(data(:,2));

%Riportiamoci al dominio [0,1]x[0,1]
data0=data;
data(:,1)=(data(:,1)-a)/(b-a);
data(:,2)=(data(:,2)-c)/(d-c);
coarse_nx=method_data.nsub_coarse(1);
coarse_ny=method_data.nsub_coarse(2);

[M,N] = size(data);
%Inizializziamo lo spazio (inizialmente è uno spazio di B-spline 
%prodotto-tensore, dunque con un solo livello)
[hmsh, hspace, geometry] = adaptivity_initialize_laplace (problem_data, method_data);
%locappr=input('Choose local approximation (1=local projector, 2=local least squares, 3=local spline/polynomial least squares):');
finest_lev=adaptivity_data.max_level;%input('Massimo numero di livelli gerarchici: ');
%fprintf('\n');

weight=ones(1,M)/ M; %inizializzazione dei pesi
tol1=adaptivity_data.tol;   %TOLLERANZA (su errore assoluto)
tol_sat=0; %flag che segnala se la tolleranza è soddisfatta 
lev=hspace.nlevels;%numero corrente di livelli dello spazio gerarchico (inizialmente =1) 
%stuck è un altro a flag per determinare se lo spazio gerarchico sia 
%invariato rispetto all'iterazione precedente
stuck=0;

%Algoritmo iterativo: ad ogni passo viene calcolata l'approssimazione ai
%minimi quadrati (con termine di penalizzazione) nello spazio gerarchico
%attuale, dopodiché si calcolano gli errori e in base ad essi si raffina la
%mesh e lo spazio gerarchico
while lev<=finest_lev && tol_sat==0 && stuck==0
    fprintf('Numero di livelli: %d\n',hspace.nlevels);    
    fprintf('Dimensione spazio: %d\n',hspace.ndof);
    clear QI_coeff
    %inizializzazione del vettore che conterrà i coefficienti della soluzione
    %QI_coeff=zeros(hspace.ndof,1); 
    disp('computing least squares solution...')  
    
    %Qui viene chiamata la funzione che calcola risolve il problema dei minimi
    %quadrati con termine di penalizzazione
    %Concretamente, vengono restituiti i coefficienti della soluzione nella base 
    %di B-spline gerarchiche
    [QI_coeff, cond_num] = getcoeff_weighted_least_squares_pen(hspace,hmsh,data,f,1e-6,weight);           

    disp('solution computed...')
    disp('evaluating solution at data points...')
    QI=sp_eval_alt (QI_coeff, hspace, data);
    
    error=abs((f-QI')); %vettore degli errori assoluti nei punti dati
    ex_tol1=find(error>tol1); %indici dei punti dove l'errore supera la tolleranza
    
    %Eventuale aggiornamento dei pesi
    weight(ex_tol1) = weight(ex_tol1).*(1+error(ex_tol1)');
    
    fprintf('Cond. number of the least squares matrix: %f\n',cond_num); %errore massimo
    fprintf('max error computed: %f\n',max(error)); %errore massimo
    fprintf('RMSE computed: %f\n\n',norm(error,2)/sqrt(length(error))); %errore quadratico medio

    % se ci sono punti in cui l'errore supera la tolleranza, si determinano
    % i supporti delle funzioni che contengono tali punti, che saranno le
    % parti della mesh da raffinare
    if numel(ex_tol1)>0
        marked=support_containing_point(hspace,hmsh,data(ex_tol1,:));
        %raffinamento della mesh e, di conseguenza, aggiornamento dello spazio
        hmsh_temp=hmsh;
        hspace_temp=hspace;
        [hmsh, hspace] = adaptivity_refine (hmsh, hspace, marked, adaptivity_data);
        lev=hspace.nlevels;
    else
        disp('tolerance satisfied')
        hmsh_temp=hmsh;
        hspace_temp=hspace;
        tol_sat=1;
    end
    if hspace.ndof==hspace_temp.ndof
        stuck=1;
    end
end

%Grafici
hmsh_plot_cells (hmsh_temp, 10, 1); %mesh
figure (1)
view(0,90)
axis off

figure(2)
scatter3(data0(:,1),data0(:,2),f); %dati originali 

figure (3) %punti originali (proiettati sul piano xy)
plot(data0(:,1),data0(:,2),'.');
axis([a b c d]) 

figure(4) %Soluzione
xplot=linspace(a,b,151);
yplot=linspace(c,d,151);
[Xplot,Yplot]=meshgrid(xplot,yplot);
QIplot=sp_eval (QI_coeff', hspace_temp, geometry,[151 151]);
QIplot=QIplot';
mesh(Xplot,Yplot,QIplot,'FaceLighting','phong','FaceColor',[.5 .5 .5],...
	'EdgeColor','none',...
	'AmbientStrength',1.0,'SpecularExponent',15,'SpecularStrength',0)   
camlight (-37.5,30+20)
axis off

figure(5) %errore
scatter3(data0(:,1),data0(:,2),error)
%% MOdel for VEgetation DYnamics in tidal salt marshes
% This is the Matlab code associate with the numerical model presented in 
% Finotello, Bertuzzo, Marani, D'Alpaos (2021):
% A minimalist model of salt-marsh vegetation dynamics driven by species competition and dispersal 
% Frontiers in Marine Science
%
%
% DISCLAIMER
% The code provided here is offered as-is, with no guarantees or technical support. 
% It is meant to provide transparency and repeatability of results and any modifications are at your own risk.
% File data for marsh elevation and vegetation distributions are also
% provided to allow user to run the model.
% The data are the same presented in the paper (see Figure 1)

close all
clear all
clc

%% COMPUTATIONAL DOMAIN SEUTP
matelev=imread('Data/Lidar_area_VoidFill.tif'); % Elevation matrix
matelev(matelev==255)=NaN;
matelev=double(matelev);
matchan=imread('Data/Mask_SanFelice_1m.tif'); % Morphological matrix
matchan=double(matchan);

[dimy, dimx]=size(matchan);

matelev=matelev.*10^-2;
matelev_ini=matelev;
matelev_noNAN=matelev;
matelev_noNAN(isnan(matelev_noNAN))=0;

% Creation of topographic domain
mattot=zeros(size(matchan));
mattot(matchan==0)=matelev(matchan==0);
mattot(matchan~=0)=NaN;

matelev=mattot;
clear mattot

%% MODEL PARAMETER
%PARAMETERS RELATED TO THE MARSH PHYSICAL PROPERTIES
pixel=1;                                    % grid resolution [m]
R=makerefmat(pixel,pixel,pixel,pixel);      % referencing matrix
zeta0=nanmean(nanmean(matelev));            % mean marsh elevation [m]
pixel2=pixel*pixel;                         % cell area [m^2]


% PARAMETERS RELATED TO VEGETATION SPECIES
N=4; % number of species
color=parula(N); % Plotting colors
color=[0,0,0;color];

% habitat quality function parameters
z0=[0.05 0.10 0.15 0.20];            % Inflection points
phi=1/(0.01*3600*24);                % fertility base rate
phi_out=phi/100;                     % fertility rate for colonization from outside the system (i.e., species that does not exist within the system at a given time)
fmax=[1 1 1 1];                      % Maximum fitness
kvege=[90 90 90 90];                 % Logistic rate

% Species competitveness
delta_max=1;
coeffang=(1-delta_max)/(abs(max(z0))-abs(min(z0)));
competitiveness=coeffang.*(z0-abs(min(z0)))+delta_max; % Species competitiveness

if isequal(repmat(N,1,4),[numel(z0),numel(fmax) numel(kvege) numel(competitiveness)])~=1
    error('Wrong number of vegetation species')
end

% Vegetation Mortality
mu0=repmat(0.005/(3600*24),1,N);       % mortality base rate [1/T]

% Dispersal Kernel for vegetation
[x,y]=meshgrid(1:dimx,1:dimy);

Disp=zeros(dimy,dimx);       % periodic boundary are used
for i=1:size(Disp,1)         % distance matrix
    for j=1:size(Disp,2)
        Disp(i,j)=sqrt((min((i-1),dimy-i+1))^2+(min((j-1),dimx-j+1))^2);
    end
end
Disp=Disp.*pixel;       %convert distance to meters

D=1;            %Kernel size
D=max(D,pixel); %WARNING: D should not be larger than pixel, otherwise numeric approximation of kernel is not satisfactory
if D/pixel>50
    D=pixel*50; %D cannot be too large otherwise numerical approximation of Kernel is not satisfactory (could be improved using non-periodic kernel but with >> computational cost)
end
Kern=exp(-Disp./D)/(2*pi*D^2)*pixel2;

if sum(sum(Kern))>=1.05 || sum(sum(Kern))<=0.95
    error('Numerical approximation error in kernel estimate')
end

Kern(1,1)=0; % to set colonization pressure on the abandoned site equal to zero (in this way, colonization pressure will depend only on the state of neighbour sites)


% % % % %%Plotting model input data
% % % % Figure
% % % % histogram(matelev)
% % % % title('Elevation distribution')
% % % % 
% % % % Figure
% % % % hold on
% % % % for j=1:N
% % % %     plot(zi,fitness(zi,z0(j),fmax(j),kvege(j)), 'Color', color(j+1,:),'LineWidth',2);
% % % % end
% % % % legend(strcat(repmat('Species ',N,1),num2str((1:N)')))
% % % % xlabel('Elevation z [m a.m.s.l.]')
% % % % ylabel('fitness f(z)')
% % % % xlim([min(matelev(:)) max(matelev(:))])
% % % % title('Vegetation Fitness')
% % % % 
% % % % Figure
% % % % hold on
% % % % for j=1:N
% % % %     plot(zi,mortality(zi,mu0(j),z0(j),fmax(j),kvege(j)),'Color', color(j+1,:),'LineWidth',1.5)
% % % % end
% % % % set(gca,'Yscale','log')
% % % % legend(strcat(repmat('Species ',N,1),num2str((1:N)'),repmat(' - z_0 = ',N,1),num2str(z0'),repmat(', f_{max} = ',N,1),num2str(fmax'),repmat(', k = ',N,1),num2str(kvege')))
% % % % xlabel('Elevation z [m a.m.s.l.]')
% % % % ylabel('Mortality rate \mu_{z}')
% % % % xlim([min(matelev(:)) max(matelev(:))])
% % % % title('Vegetation Mortality')


%% SIMULATION PARAMETERS
Ttot=50*365*24*3600;     % Total duration of eco-morphodynamic simulation [seconds]
dt_sim=(365/12)*24*3600; % Time-step for eco-morphodynamic simulations [seconds]
dt_print=365*24*3600;    % Time-step for printing results [seconds]

HabitatQuality_mat=zeros(dimy,dimx,N);
for j=1:N
    HabitatQuality_mat(:,:,j)=fitness(matelev_noNAN,z0(j),fmax(j),kvege(j));
end

sim_repeat=20;  %number of simulation repetitions (to filter out fluctuations due to demographic stochasticity)
sim_count=0;    %simulation counting

%% RUN SIMULATION
while sim_count<sim_repeat
    sim_count=sim_count+1;
    disp(strcat('Simulation n.',num2str(sim_count),'\', num2str(sim_repeat)))
    
    % initial vegetation state
    mat_status=randi([0 N], size(matelev)); %random vegetatio distribution
    mat_status(isnan(matelev))=NaN;
    
    monitor=zeros(floor(Ttot/dt_print)+1,1);        %to monitor sim results
    species_perc=zeros(floor(Ttot/dt_print)+1,N+1); %to monitor sim results

    for t=0:dt_sim:Ttot
%         disp(num2str(t/(365*24*3600)))
        
        if mod(t,dt_print)==0
            for j=0:N
                species_perc(floor(t/dt_print)+1,j+1)=sum(mat_status(:)==j)/(sum(~isnan(mat_status(:))))*100;
            end
        end
        
        %Colonization rate [1/s]
        C=zeros(dimy,dimx,N);
        for j=1:N
            C(:,:,j)=competitiveness(j)*HabitatQuality_mat(:,:,j).*(ifft2(fft2(Kern).*fft2((mat_status==j).*squeeze(HabitatQuality_mat(:,:,j))*phi))+phi_out);
        end
        
        %Vegetation mortality
        mat_mortality=zeros(size(mat_status));
        for j=1:N
            iiddxx=mat_status==j;
            mat_mortality(iiddxx)=mortality(matelev(iiddxx),mu0(j),z0(j),fmax(j),kvege(j)); %morality rate
        end
        mat_status(rand(size(mat_mortality))<1-exp(-mat_mortality*dt_sim))=0;
        
        %Vegetation growback
        idx=find(mat_status==0); % non vegetated pixels
        
        %compute total colonizaton rate
        cumulative_colon_rate=cumsum(C,3);
        total_colon_rate=cumulative_colon_rate(:,:,end);
        i_p_colon=idx(rand(size(idx))<1-exp(-total_colon_rate(idx)*dt_sim));  %index of  pixels that will be colonized
        
                % select which species colonize the i_p_colon pixels
        prob_mat=cumulative_colon_rate./total_colon_rate;
        
        selection=rand(length(i_p_colon),1);
        mat_status(i_p_colon)=N;
        for j=N-1:-1:1
            col=selection<prob_mat((j-1)*(dimx*dimy)+i_p_colon);
            mat_status(i_p_colon(col))=j;
        end
    end
    
    %Save species percent for sim_count
    save(strcat('Species_Percent_', num2str(sim_count),'.mat'), 'species_perc')
    save(strcat('mat_status',num2str(sim_count),'.mat'), 'mat_status')
end

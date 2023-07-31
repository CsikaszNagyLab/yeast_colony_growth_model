function [p, p3, popnum, final_population, dead_cells, area_matrix]=yeast_model(inputname, outputname)
%yeast growth model 

tstart=tic; %measure time

%read in general parameters:
dt=dlmread(inputname, ';', 'B1..B1')
dx=dlmread(inputname, ';', 'B2..B2')
Dp=dlmread(inputname, ';', 'B3..B3') %top nutrient layer convolution coefficient
Dp2=dlmread(inputname, ';', 'B4..B4') %deep nutrient layer convolution coefficient
Dp3=dlmread(inputname, ';', 'B5..B5') %ECM or signal convolution coefficient
mx_decay=dlmread(inputname, ';', 'B6..B6'); % ECM or signal decay parameter
max_x=dlmread(inputname, ';', 'B7..B7')
max_y=dlmread(inputname, ';', 'B8..B8')
boundary_cond=dlmread(inputname, ';', 'B9..B9')
nt=dlmread(inputname, ';', 'B10..B10')
recharge_rate=dlmread(inputname, ';', 'B11..B11')
init_nutrient=dlmread(inputname, ';', 'B12..B12')
init_nutrient2=dlmread(inputname, ';', 'B13..B13')

%cell parameters: (max 24 species)
init_E=dlmread(inputname, ';', 'B16..Y16') 
p_uptake=dlmread(inputname, ';', 'B17..Y17')
p_uptake_eff=dlmread(inputname, ';', 'B18..Y18')
div_th=dlmread(inputname, ';', 'B19..Y19')
div_nbor=dlmread(inputname, ';', 'B20..Y20')
slope1=dlmread(inputname, ';', 'B21..Y21') %division distance decrease: give high value for no decrease
death_th=dlmread(inputname, ';', 'B22..Y22')
metab_E=dlmread(inputname, ';', 'B23..Y23')
colony_type=dlmread(inputname, ';', 'B24..Y24')
g0_factor=dlmread(inputname, ';', 'B25..Y25')  %G0 state nutrient_uptake
g0_th=dlmread(inputname, ';', 'B26..Y26')
extracell_mx_prod=dlmread(inputname, ';', 'B27..Y27') %extracell matrix production
exrtracell_G0=dlmread(inputname, ';', 'B28..Y28') %G0 extracell matrix production
a=dlmread(inputname, ';', 'B29..Y29') %signal effect parameter (rate change)
c=dlmread(inputname, ';', 'B30..Y30') %signal effect parameter (x=c when f(x)==0.5)
branchprob=dlmread(inputname, ';', 'B31..Y31')
div_distrdev=dlmread(inputname, ';', 'B32..Y32')
div_distrdev(find(div_distrdev))=pi./div_distrdev(find(div_distrdev));
initnum_species=dlmread(inputname, ';', 'B33..Y33')
initdev=dlmread(inputname, ';', 'B34..Y34')
epow=dlmread(inputname, ';', 'B35..Y35') %cell distribution inside the initial area (circle) of cells: epow=0.25 -> more cells around the edge, 0.5 -> uniform distribution, 1 -> more cells around the centre 
initpos(1,:)=dlmread(inputname, ';', 'B36..Y36')
initpos(2,:)=dlmread(inputname, ';', 'B37..Y37')
initpos=initpos';

%visualisation parameters:
isdraw=dlmread(inputname, ';', 'B40..B40')
plotstep=dlmread(inputname, ';', 'B41..B41')


%initialise:
x = 1:max_x;
y = 1:max_y;
px = length(x);
py = length(y);
if init_nutrient==0
init_nutrient=1e-100;
end

sp_id=find(initnum_species);
species_number=sp_id(end);
species_name=1:species_number

area_matrix = zeros(nt+1,species_number); %area measurement
popnum=zeros(nt+1,species_number+1+species_number*2); %store all the states (active, G0)
popnum(1,:)=[sum(initnum_species), initnum_species(1:species_number), initnum_species(1:species_number), zeros(1,species_number)]; %initial states
species_pos=nan(popnum(1,1),2);
species_E=nan(popnum(1,1),1);
species=nan(popnum(1,1),1);
species_div=nan(popnum(1,1),1);
species_divpref=nan(popnum(1,1),1);
borders=[0 cumsum(initnum_species)];

p = init_nutrient*ones(px,py);%top nutrient layer
p2 = init_nutrient2*ones(px,py);%deep nutrient layer
p3=zeros(px,py);%extracellular matrix (ECM) os signal layer

%inhomogeneous initial nutrient distribution
%{
p(px/5+1:end,:)=zeros(px*4/5,py);
p2(px/5+1:end,:)=zeros(px*4/5,py);
%}

%initalise agents:
for i=1:species_number
theta = rand(initnum_species(i),1)*(2*pi);
r = power(rand(initnum_species(i),1),epow(i))*initdev(i);
species_pos(borders(i)+1:borders(i+1),:)=[initpos(i,1) + r.*cos(theta), initpos(i,2) + r.*sin(theta)];

    xn=find(species_pos(borders(i)+1:borders(i+1),1)>max_x);
    xk=find(species_pos(borders(i)+1:borders(i+1),1)<0);
    yn=find(species_pos(borders(i)+1:borders(i+1),2)>max_y);
    yk=find(species_pos(borders(i)+1:borders(i+1),2)<0);
    if boundary_cond==0
    %periodic boundary:
    species_pos(xn,1)=species_pos(xn,1)-max_x;
    species_pos(xk,1)=max_x+species_pos(xk,1);
    species_pos(yn,2)=species_pos(yn,2)-max_y;
    species_pos(yk,2)=max_y+species_pos(yk,2);
    elseif boundary_cond==1
    %fix boundary: 
    species_pos(xn,1)=max_x;
    species_pos(xk,1)=0.01;
    species_pos(yn,2)=max_y;
    species_pos(yk,2)=0.01;    
    end
    
species(borders(i)+1:borders(i+1))=species_name(i)*ones(initnum_species(i),1);
%species_E(borders(i)+1:borders(i+1))=init_E(i).*ones(initnum_species(i),1); %same initial energy for all cells -> synchronised cell poulation
%species_E(borders(i)+1:borders(i+1))=-0.5.*log(exp(-init_E(i)*2)+(exp(-div_th(i)*2)-exp(-init_E(i)*2)).*rand([initnum_species(i),1])); %exponential initial energy distribution
species_E(borders(i)+1:borders(i+1))=init_E(i)+((div_th(i)-init_E(i)).*rand(initnum_species(i),1)); %uniform inital energy distribution
species_div(borders(i)+1:borders(i+1))=ones(initnum_species(i),1);

if colony_type(i)==1
species_divpref(borders(i)+1:borders(i+1))=2*pi*rand(initnum_species(i),1);
end
end

species_grid=ceil(species_pos);

%area calculation: (number of diffusion cells with agents)
for i=1:species_number
    area_matrix(1,i)=size(unique(species_grid(species==species_name(i),:), 'row'),1);
end

dead_species_pos=[];
dead_species_grid=[];
dead_species_step=[];
dead_species=[];

%plot initial state:
plotcolors={rgb('LightGreen'), rgb('Green'), rgb('LightSalmon'), rgb('Red'), rgb('LightBlue'), rgb('Blue'), rgb('Lavender'), rgb('Purple'),rgb('LightGray'), rgb('Gray'), rgb('Cornsilk'), rgb('Brown')};
%plotcolors={'g.', 'g.', 'r.', 'r.', 'c.', 'o.', 'y.', 'b.'};
%plotcolors={'k.', 'g.', 'w.', 'y.', 'c.', 'o.', 'r.', 'b.'};
%plotcolors={'r.', 'g.'};

if isdraw==1
f1=figure;
    %subplot(2,1,2)
    %imagesc(p3);
    %axis image;
    %subplot(2,1,1)
    imagesc(p+p2);
    caxis([0 init_nutrient+init_nutrient2]); 
    axis image;
hold on;
for i=1:species_number
    %plot G0 with different colors:
        %plot(species_pos((species==species_name(i) & species_div==0),2),species_pos((species==species_name(i) & species_div==0),1), '.', 'color', plotcolors{(i*2)-1}, 'MarkerSize', 1);
        %plot(species_pos((species==species_name(i) & species_div==1),2),species_pos((species==species_name(i) & species_div==1),1), '.', 'color', plotcolors{(i*2)}, 'MarkerSize', 1);
    plot(species_pos(species==species_name(i),2),species_pos(species==species_name(i),1), 'r.', 'MarkerSize', 1);
end
hold off;
title('plate diffusion time: 0');
waitforbuttonpress; %start the simulation for button press
end


%simulation steps:
for t=1:dt:nt 
    
    %{
    %wet-dry change:
    if t==169 %change time
        Dp_ratio=Dp2/Dp;
        Dp=30
        Dp2=Dp*Dp_ratio
        div_nbor(1)=0.63
    end
    %}
    
    %life cycle of the cells:
    new_cell=zeros(popnum(t,1),1);
    survive_cell=zeros(popnum(t,1),1);
    new_dir=nan(popnum(t,1),1);
    
    order=randperm(popnum(t,1)); %reorder agents for avoiding systematic bias in nutrient availability 
    
    species_pos=species_pos(order,:);
    species_grid=species_grid(order,:);
    species_E=species_E(order,:);
    species=species(order,:);
    species_div=species_div(order,:);
    species_divpref=species_divpref(order,:);
    
    %agent's life circle:
    for s=1:popnum(t,1) 
        s_akt=species(s,1); %ID of the current agent  
    
    %nutrient uptake:
    %if there is enough nutrient (at the current iteration in the grid position of the cell), the cell takes it, if not, it consumes the leftover from the upper nutrent layer 
    %metabolism:
    if species_div(s)==0 %G0 metabolism: less nutrient uptake and less energy consumption 
        p_uptake_akt=p_uptake(s_akt)*g0_factor(s_akt);
        metab_e_akt=metab_E(s_akt)*g0_factor(s_akt);
        extracell_mx_akt=extracell_mx_prod(s_akt)*exrtracell_G0(s_akt); %ECM or signal production
    elseif species_div(s)==1 %normal cell circle: normal metabolism 
        p_uptake_akt=p_uptake(s_akt);
        metab_e_akt=metab_E(s_akt);
        extracell_mx_akt=extracell_mx_prod(s_akt); %ECM or signal production
    end
   
    %if there is not enough nutrient:
    if p(species_grid(s,1), species_grid(s,2))<p_uptake_akt
        p_uptake_akt=p(species_grid(s,1), species_grid(s,2));
    end
    
    %update nutrient, ECM/signal layers and agent energy level
    p(species_grid(s,1), species_grid(s,2))=p(species_grid(s,1), species_grid(s,2))-p_uptake_akt;
    species_E(s)=species_E(s)+p_uptake_akt*p_uptake_eff(s_akt)-metab_e_akt;
    p3(species_grid(s,1), species_grid(s,2))=p3(species_grid(s,1), species_grid(s,2))+extracell_mx_akt;
        
   %decide state: G0 or active (division):
    %high signal concentration switch the cells to G0 state (with sigmoid characteristics), additional constraint 
    %y = (b+d)./(1 + exp(-a.*(x-c)))+d; 
    %parameters: 'c' corresponds to the x value where y = 0.5, x axis shift ;'a' is the rate of change ;'b' is the max value, high of the plot; 'd' is y axis shift
   if species_div(s)==1 && (species_E(s)<g0_th(s_akt) || rand<1/(1+exp(-a(s_akt).*(p3(species_grid(s,1), species_grid(s,2))-c(s_akt)))))
       species_div(s)=0; %G0
       %species_E(s)=death_th(s_akt); %toxin kills the cells (it make the simulation faster) 
   elseif colony_type(s_akt)==0 && species_div(s)==0 && species_E(s)>=g0_th(s_akt) && rand>=1/(1+exp(-a(s_akt).*(p3(species_grid(s,1), species_grid(s,2))-c(s_akt))))
       species_div(s)=1;
   end
        
    %cell division:
    if colony_type(s_akt)==1 %filamentous growth
        if species_E(s)>=div_th(s_akt)
            div=rand;
            if species_div(s)==0 && div<branchprob(s_akt)
                species_E(s)=species_E(s)-init_E(s_akt); %asymmetric cell division  
                new_cell(s,1)=s_akt;
                new_dir(s,1)=2*pi*rand;
            elseif species_div(s)==1 && div>0
                species_E(s)=species_E(s)-init_E(s_akt); %asymmetric cell division
                new_cell(s,1)=s_akt;
                new_dir(s,1)=species_divpref(s)+div_distrdev(s_akt).*randn;
                species_div(s)=0; 
            end        
        end 
    else %colony forming growth
        if species_E(s)>=div_th(s_akt) && species_div(s)==1
            species_E(s)=species_E(s)-init_E(s_akt); %asymmetric cell division
            new_cell(s,1)=s_akt;
            new_dir(s,1)=2*pi*rand;
        end
    end
    
    %cell death:
    %
    if species_E(s)>death_th(s_akt) 
       survive_cell(s,1)=1; %survivers
    else %not survive
       p(species_grid(s,1), species_grid(s,2))=p(species_grid(s,1), species_grid(s,2))+species_E(s)*1/p_uptake_eff(s_akt); %(ration of the) remaining nutrient of the cell flows back to the nutrient layer
    end
    %}  
    %if the signal kills the cells
    %{
    if species_E(s)>death_th(s_akt) && rand>1/(1+exp(-a(s_akt).*(p3(species_grid(s,1), species_grid(s,2))-c(s_akt))))
       survive_cell(s,1)=1; 
    else
       p(species_grid(s,1), species_grid(s,2))=p(species_grid(s,1), species_grid(s,2))+species_E(s)*1/p_uptake_eff(s_akt); %(ratio of the) remaining nutrient of the cell flows back to the nutrient layer
    end
    %}
    
    end
 
    %new cells:
    if new_cell==0
    new_id=double.empty(0,1);
    else
    new_id=find(new_cell);
    end
    new_species=new_cell(new_id);
    new_dir=new_dir(new_id);
    new_pos=[species_pos(new_id,1)+div_nbor(new_species)'.*cos(new_dir).*(1-t./(t+slope1(s_akt))), species_pos(new_id,2)+div_nbor(new_species)'.*sin(new_dir).*(1-t./(t+slope1(s_akt)))]; %time dependent cell division, older colonies can dry up -> shorter division distances: 1-t./(t+slope1); 
    xn=find(new_pos(:,1)>max_x);
    xk=find(new_pos(:,1)<0);
    yn=find(new_pos(:,2)>max_y);
    yk=find(new_pos(:,2)<0);
    if boundary_cond==0
    %periodikus határ
    new_pos(xn,1)=new_pos(xn,1)-max_x;
    new_pos(xk,1)=max_x+new_pos(xk,1);
    new_pos(yn,2)=new_pos(yn,2)-max_y;
    new_pos(yk,2)=max_y+new_pos(yk,2);
    elseif boundary_cond==1
    %fix boundary
    new_pos(xn,1)=max_x;
    new_pos(xk,1)=0.01;
    new_pos(yn,2)=max_y;
    new_pos(yk,2)=0.01;
    end
    
    %update agent list: delete dead cells
    dead_id=find(~survive_cell);
    dead_species_pos=[dead_species_pos; species_pos([dead_id],:)];
    dead_species_grid=[dead_species_grid; species_grid([dead_id],:)];
    dead_species=[dead_species; species(dead_id)];
    dead_species_step=[dead_species_step; t.*ones(size(dead_id,1),1)];
    
    survive_id=find(survive_cell);
    species_pos=species_pos([survive_id],:);
    species_grid=species_grid((survive_id),:);
    species_E=species_E((survive_id),:);
    species=species((survive_id),:);
    species_div=species_div((survive_id),:);
    species_divpref=species_divpref((survive_id),:);
    
    %update agent list: add new daughter cells
    species_pos=[species_pos; new_pos];
    species_grid=[species_grid; ceil(new_pos)];
    species_E=[species_E; init_E(1).*ones(size(new_id,1),1)];  
    species=[species; new_species];
    species_div=[species_div; ones(size(new_id,1),1)];
    species_divpref=[species_divpref; new_dir];
    
    %{
    %put new cells no the plate:
    if t==480
        i=1;
        pos=[310, 400];
        theta = rand(initnum_species(i),1)*(2*pi);
        r = power(rand(initnum_species(i),1),epow(i))*initdev(i);
        new_pos2=[pos(1) + r.*cos(theta), pos(2) + r.*sin(theta)];
        species_pos=[species_pos; new_pos2];
        species_grid=[species_grid; ceil(new_pos2)];
        species_E=[species_E; init_E(1).*ones(size(new_pos2,1),1)];  
        species=[species; species_name(1).*ones(size(new_pos2,1),1)];
        species_div=[species_div; ones(size(new_pos2,1),1)];
        species_divpref=[species_divpref; 2*pi*rand(initnum_species(i),1)];
    end
    %}
    
    sp_states=[sum(species==species_name), sum(species==species_name & species_div==1), sum(species==species_name & species_div==0)]; %iterate through strains sum up and store the numbers of active (1) and G0 (0) cells 
    popnum(t+1,:)=[size(species_E,1), sp_states]; %actual cell and cell state numbers for each strain
    
    %area calculation (diffusion cells occupied by agents):
    for i=1:species_number
        %area_matrix(t+1,i)=size(unique(species_grid(species==species_name(i),:), 'row'),1); %only living cells   
        area_matrix(t+1,i)=size(unique([species_grid(species==species_name(i),:); dead_species_grid(dead_species==species_name(i),:)], 'row'),1); %area with dead cells (all area ever)
    end
    
    
    %material distribution with convolution: (mimics diffusion)
    p=imgaussfilt(p, Dp); % Gaussian filtering for nutrient distribution
    p2=imgaussfilt(p2, Dp2); % Gaussian filtering for nutrient distribution
    p3=imgaussfilt(p3, Dp3); % Gaussian filtering for ECM/signal distribution
    p3=(1-mx_decay)*p3; %ECN/signal decay

    %material transfer between layers:
    recharge=(p2/(init_nutrient2/init_nutrient)-p)*recharge_rate;
    p=p+recharge;
    p2=p2-recharge;
    
    %{
    % nutrient supply at given steps: (p2: lower/deep layer)
    supplystep=168;
    supply_amount=7;
    if mod(t, supplystep)==0
    p2(1,:) = p2(1,:)+supply_amount*init_nutrient2*ones(1,py);
    p2(end,:) = p2(end,:)+supply_amount*init_nutrient2*ones(1,py);
    p2(2:end-1,1) = p2(2:end-1,1)+supply_amount*init_nutrient2*ones(px-2,1);
    p2(2:end-1,end) = p2(2:end-1,end)+supply_amount*init_nutrient2*ones(px-2,1);
    end
    %}
    
    %plot the actual state
    if mod(t,plotstep)==0 
    [t, popnum(t+1, :), area_matrix(t+1, :)] %print out the strain and states numbers and the areas
    if isdraw==1
    pause(0.001);
  
    %subplot(2,1,2)
    %imagesc(p3);
    %axis image;
    %colorbar; 
    %subplot(2,1,1)
    imagesc(p+p2);
    caxis([0 init_nutrient+init_nutrient2]);
    axis image;
    hold on;
    for i=1:species_number
    %plot G0 with different colors:
    %plot(species_pos((species==species_name(i) & species_div==0),2),species_pos((species==species_name(i) & species_div==0),1), '.', 'color', plotcolors{(i*2)-1}, 'MarkerSize', 1);
    %plot(species_pos((species==species_name(i) & species_div==1),2),species_pos((species==species_name(i) & species_div==1),1), '.', 'color',  plotcolors{(i*2)}, 'MarkerSize', 1);
    plot(species_pos(species==species_name(i),2),species_pos(species==species_name(i),1), 'r.', 'MarkerSize', 1);
    end
    grid on;
    hold off;
    title(['plate diffusion time: ', num2str(t)]);
        %waitforbuttonpress;
    end
    end
       
    %{
    %save at given time points
    if any(t==[48, 312])
    %plot cell density: number of cells/grid positions in a given timestep
    a1=accumarray([[species_grid(:,1); dead_species_grid(:,1)],[species_grid(:,2); dead_species_grid(:,2)]],1);
    a2=zeros(size(p));
    a2(1:size(a1,1),1:size(a1,2))=a1;
    %a2=imgaussfilt(a2,2);
    f2=figure;
    %surf(a2);
    %shading interp;
    imagesc(a2);
    colorbar;
    axis image
    print([outputname, '_t', num2str(t), '_colony_allcells'],'-djpeg'); %'_neardot_colony_allcells'
    save([outputname, '_t', num2str(t), '.mat'], 'p', 'popnum', 'species', 'species_pos', 'species_grid', 'species_E', 'species_div', 'species_divpref', 'dead_species', 'dead_species_pos', 'dead_species_grid', 'dead_species_step', 'area_matrix');

    figure(f1);
    end
    %}
    
end


%generating and saving the output:
final_population=[species, species_pos, species_grid, species_E, species_div, species_divpref];
dead_cells=[dead_species, dead_species_pos, dead_species_grid, dead_species_step];  
if isdraw==1
print([outputname, '_2colony_final'],'-djpeg');
end
%
%save: in .mat
save([outputname, '.mat'], 'p', 'popnum', 'species', 'species_pos', 'species_grid', 'species_E', 'species_div', 'species_divpref', 'dead_species', 'dead_species_pos', 'dead_species_grid', 'dead_species_step', 'area_matrix');

%cell states, occupied areas for each species in each time steps: 
sn=textscan(num2str(species_name),'%s');
my_header1=horzcat({'step', 'sum'}, sn{:}');
fmt=repmat('%s\t',1,length(my_header1));
fmt(end)='n';
fid=fopen([outputname, '_popnum.tsv'],'wt');
fprintf(fid, fmt, my_header1{:});
fclose(fid);
dlmwrite([outputname, '_popnum.tsv'], [(0:nt)', popnum, area_matrix], '-append', 'delimiter', '\t', 'newline', 'pc');

%final population:
my_header2={'species', 'position x', 'position y', 'grid x', 'grid y', 'energy', 'divition state', 'division preferency'};
fmt=repmat('%s\t',1,length(my_header2));
fmt(end)='n';
fid=fopen([outputname, '_final_cell_properties.tsv'],'wt');
fprintf(fid, fmt, my_header2{:});
fclose(fid);
dlmwrite([outputname, '_final_cell_properties.tsv'], final_population, '-append', 'delimiter', '\t', 'newline', 'pc');

%dead cells:
my_header3={'species', 'position x', 'position y', 'grid x', 'grid y', 'death step'};
fmt=repmat('%s\t',1,length(my_header3));
fmt(end)='n';
fid=fopen([outputname, '_dead_cells.tsv'],'wt');
fprintf(fid, fmt, my_header3{:});
fclose(fid);
dlmwrite([outputname, '_dead_cells.tsv'], dead_cells, '-append', 'delimiter', '\t', 'newline', 'pc');
%}

elapsedtime=toc(tstart)/60
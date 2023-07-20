%run yeast growth model:
inputname='1species_parameters.csv';%input parameter file 
outputname='output/experiment1';%folder and name of the experiment 
[p, p3, popnum, final_population, dead_cells, area_matrix]=yeast_model(inputname, outputname); %rum the model

%output visualizations:

%%
%cell numbers:
figure;
nt=size(popnum,1)-1;
hold on;
plot(0:nt, popnum(:,2), '-', 'color', rgb('Green')); 
if size(popnum,2)>4
plot(0:nt, popnum(:,3), '-', 'color', rgb('Red')); 
end
xlabel('steps');
ylabel('population number');
legend({'strain 1', 'strain 2'});
print([outputname, '_popnum'],'-djpeg');

%%
% colony area (occupied cells):
figure;
nt=size(popnum,1)-1;
hold on;
plot(0:nt, area_matrix(:,1), '-', 'color', rgb('Green')); 
if size(area_matrix,2)>1
plot(0:nt, area_matrix(:,2), '-', 'color', rgb('Red')); 
end
legend({'species 1', 'species 2'},'Location','northwest');
xlabel('steps');
ylabel('colony area');
print([outputname, '_colony_area'],'-djpeg');

%%
%plot cell density: cells/grid position

%living plus dead cells:
a=accumarray([[final_population(:,4); dead_cells(:,4)],[final_population(:,5); dead_cells(:,5)]],1);
a2=zeros(size(p));
a2(1:size(a,1),1:size(a,2))=a;
%a2=imgaussfilt(a2,2);
figure;
%surf(a2);
%shading interp;
imagesc(a2);
colorbar;
axis image
print([outputname, '_colony_allcells'],'-djpeg');

%living (active plus G0) cells:
a3=accumarray([final_population(:,4),final_population(:,5)],1);
a4=zeros(size(p));
a4(1:size(a3,1),1:size(a3,2))=a3;
%a4=imgaussfilt(a4,1);
figure;
%surf(a4);
%shading interp;
imagesc(a4);
colorbar;
axis image;
print([outputname, '_colony'],'-djpeg');

%active cells:
a3=accumarray([final_population(final_population(:,7)==1,4),final_population(final_population(:,7)==1,5)],1);
a4=zeros(size(p));
a4(1:size(a3,1),1:size(a3,2))=a3;
%a4=imgaussfilt(a4,1);
figure;
%surf(a4);
%shading interp;
imagesc(a4);
colorbar;
axis image;
caxis([0 30]);
print([outputname, '_colony_active'],'-djpeg');

%%
%{
%show only one strain (e.g. strain 1):

%living plus dead cells:
species=1:(size(popnum,2)-2);
a=accumarray([[final_population(final_population(:,1)==1,4); dead_cells(dead_cells(:,1)==1,4)],[final_population(final_population(:,1)==1,5); dead_cells(dead_cells(:,1)==1,5)]],1);
a2=zeros(size(p));
a2(1:size(a,1),1:size(a,2))=a;
    %a2=imgaussfilt(a2,4);
figure;
%surf(a2);
%shading interp;
imagesc(a2);
colorbar;
print([outputname, '_colony_allcells_sp1'],'-djpeg');
size(find(a2),1)

%living (active plus G0) cells:
a=accumarray([final_population(final_population(:,1)==1,4),final_population(final_population(:,1)==1,5)],1);
a2=zeros(size(p));
a2(1:size(a,1),1:size(a,2))=a;
    %a2=imgaussfilt(a2,4);
figure;
%surf(a2);
%shading interp;
imagesc(a2);
colorbar;
print([outputname, '_colony_sp1'],'-djpeg');

%active cells:
a=accumarray([final_population((final_population(:,1)==1 & final_population(:,7)==1),4),final_population((final_population(:,1)==1 & final_population(:,7)==1),5)],1);
a2=zeros(size(p));
a2(1:size(a,1),1:size(a,2))=a;
    %a2=imgaussfilt(a2,4);
figure;
%surf(a2);
%shading interp;
imagesc(a2);
colorbar;
print([outputname, '_colony_active_sp1'],'-djpeg');
%}

%%
%{
%final population plot:
figure;
imagesc(p);
    caxis([0 100]);
    hold on;
    %plot(species_pos(species==1,2),species_pos(species==1,1), 'g.', 'MarkerSize', 1);
    %plot(species_pos(species==2,2),species_pos(species==2,1), 'r.', 'MarkerSize', 1); 
    hold off;
%}

%%
%
%show final populations blured:
load([outputname, '.mat']); %load results from output folder
sigmagauss=2; %bluring parameter
figure;
%a=accumarray([species_grid(species(:,1)==1,1),species_grid(species(:,1)==1,2)],1); %only living strain 1 cells  
a=accumarray([species_grid(:,1), species_grid(:,2); dead_species_grid(:,1), dead_species_grid(:,2)],1); %with dead cells

a2=zeros(size(p));
a2(1:size(a,1),1:size(a,2))=a;
    a2=imgaussfilt(a2,sigmagauss);
    a2=a2/max(max(a2));
my_image=zeros(size(p,1),size(p,2),3);
my_image(:,:,1)=a2;
my_image(:,:,2)=a2;
my_image(:,:,3)=a2;
imshow(my_image);
print([outputname, '_image_colony_sp1'],'-djpeg');
%}

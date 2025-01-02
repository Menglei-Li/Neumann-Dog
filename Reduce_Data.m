% Handle FEM data and FFT data

clear; clc;

% RVE model resolution and size
RVE_N1 = 1000; RVE_N2 = 120; RVE_N3 = 2;
RVE_L1 = 200; RVE_L2 = 24.00; RVE_L3 = 2.0;

zero_matrix = zeros(RVE_N1,RVE_N2,RVE_N3);
unit_matrix = ones(RVE_N1,RVE_N2,RVE_N3);

RVE_Label = Generate_Marco_non_periodic_Dog(RVE_L1,RVE_L2,RVE_N1,RVE_N2,RVE_N3);
RVE_Label(:,1:6,:) = 2; RVE_Label(:,7:10,:) = 0;
RVE_Label(:,111:120,:) = 2; RVE_Label(:,111:114,:) = 0;
RVE_Label(1:50,:,:) = 2;
RVE_Label(RVE_N1-49:RVE_N1,:,:) = 2;


%%%% Extract the two arrays (element number and field) extracted from the ODB file and re-characterize them as a two-dimensional array and visualize them %%%%
dataS11 = load('FEM\S11.mat').data;
dataS22 = load('FEM\S22.mat').data;
dataS12 = load('FEM\S12.mat').data;
dataLabel = load('FEM\Label.mat').data;

% If the FEM model uses reduced integration elements, no averaging is required (one integration point per element)

% --------------------- The order of unit numbering in FEM model is X first and then Y ---------------------
RVE_N1 = 900; RVE_N2 = 120; RVE_N3 = 2;

FEM_sig11 = zeros(RVE_N1,RVE_N2,RVE_N3);
FEM_sig22 = zeros(RVE_N1,RVE_N2,RVE_N3);
FEM_sig12 = zeros(RVE_N1,RVE_N2,RVE_N3);

for n_data = 1:size(dataS11,2)
    n_x = floor( double(dataLabel(n_data)) / double(RVE_N2*RVE_N3));
    nn_x = n_x;
    if dataLabel(n_data)-nn_x*(RVE_N2*RVE_N3)==0
        n_x = n_x-1;
        nn_x = n_x;
    end
    n_x = n_x+1;

    n_z = floor( double(dataLabel(n_data)-nn_x*(RVE_N2*RVE_N3)) / double(RVE_N2) );
    nn_z = n_z;
    if dataLabel(n_data)-nn_x*(RVE_N2*RVE_N3)-nn_z*RVE_N2==0
        n_z = n_z-1;
        nn_z = n_z;
    end
    n_z = n_z+1;

    n_y = dataLabel(n_data)-nn_x*(RVE_N2*RVE_N3)-nn_z*RVE_N2;
    if dataLabel(n_data)-nn_x*(RVE_N2*RVE_N3)-nn_z*RVE_N2==0
        n_y = RVE_N2;
    end

    FEM_sig11(n_x,n_y,n_z) = dataS11(n_data);
    FEM_sig22(n_x,n_y,n_z) = dataS22(n_data);
    FEM_sig12(n_x,n_y,n_z) = dataS12(n_data);

end

% Extract FFT calculation results
FFT_sig11 = load('FFT_sig11.mat').RVE_sig11; FFT_sig11 = FFT_sig11(51:RVE_N1+50,:,:);
FFT_sig22 = load('FFT_sig22.mat').RVE_sig22; FFT_sig22 = FFT_sig22(51:RVE_N1+50,:,:);
FFT_sig12 = load('FFT_sig12.mat').RVE_sig12; FFT_sig12 = FFT_sig12(51:RVE_N1+50,:,:);

subplot(1, 3, 1);
imagesc(squeeze(FEM_sig11(:,:,1)));
colormap(jet);
colorbar('northoutside');
set(gca, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none');

subplot(1, 3, 2);
imagesc(squeeze(FFT_sig11(:,:,1)));
colormap(jet);
colorbar('northoutside');
set(gca, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none');

subplot(1, 3, 3);
imagesc( squeeze( (FEM_sig11(:,:,1)-FFT_sig11(:,:,1)) ) );
colormap(jet);
colorbar('northoutside');
set(gca, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none');

DIF_sig11 = FEM_sig11-FFT_sig11; RVE_Label = RVE_Label(51:RVE_N1+50,:,:);
DIF_sig11(RVE_Label==0)=-100;
% Output all results as .Vtk files
% Set the working directory
RVE_folder_load = strcat('FFT Dog');
if ~exist(RVE_folder_load,'dir')
    mkdir(RVE_folder_load); 
end
RVE_file_path_load = strcat(pwd,'\',RVE_folder_load,'\');

Marco_Vtk(FEM_sig11,RVE_file_path_load,'FEM_sig11','S11',1);
Marco_Vtk(FFT_sig11,RVE_file_path_load,'FFT_sig11','S11',1);
Marco_Vtk(DIF_sig11,RVE_file_path_load,'DIF_sig11','S11',1);

% Extract the results for 5 positions (P1-P5)

% P1
FFT_P1 = FFT_sig11(1,11:110,1)'; FEM_P1 = FEM_sig11(1,11:110,1)';

% P2
FFT_P2 = FFT_sig11(164,11:110,1)'; FEM_P2 = FEM_sig11(164,11:110,1)';

% P3
FFT_P3 = FFT_sig11(225,11:110,1)'; FEM_P3 = FEM_sig11(225,11:110,1)';

% P4
FFT_P4 = FFT_sig11(292,11:110,1)'; FEM_P4 = FEM_sig11(292,11:110,1)';

% P5
FFT_P5 = FFT_sig11(450,11:110,1)'; FEM_P5 = FEM_sig11(450,11:110,1)';

save("P1_P5.mat","FFT_P1","FEM_P1","FFT_P2","FEM_P2","FFT_P3","FEM_P3","FFT_P4","FEM_P4","FFT_P5","FEM_P5");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% functions used %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Macro_Label = Generate_Marco_non_periodic_Dog(Marco_L1,Marco_L2,Marco_N1,Marco_N2,Marco_N3)

    Macro_Label=zeros(Marco_N1,Marco_N2,Marco_N3);
    d_x=Marco_L1/(Marco_N1);d_y=Marco_L2/(Marco_N2);

    % 规定中心板子
    c_x1 = 62.5; c_y1 = -85.5;
    c_x2 = 62.5; c_y2 = 109.5;
    c_x3 = 117.5; c_y3 = -85.5;
    c_x4 = 117.5; c_y4 = 109.5; r = 92.5;

    Macro_Label(1:50,:,:) = 2;
    Macro_Label(Marco_N1-49:Marco_N1,:,:) = 2;

    for xi=51:Marco_N1-50
        n_x=(xi-50)*d_x-0.5*d_x;
        for yi=1:Marco_N2
            n_y=yi*d_y-0.5*d_y;
            if n_y>=7 && n_y<=17
                Macro_Label(xi,yi,:)=1;
            elseif n_x<=32.5 || n_x>=147.5
                if n_y>=2 && n_y<=22
                    Macro_Label(xi,yi,:)=1;
                end
            elseif n_x>32.5 && n_x<62.5 && n_y>=2 && n_y<=7
                dis1 = sqrt((n_x-c_x1)^2+(n_y-c_y1)^2);
                if dis1>=r
                    Macro_Label(xi,yi,:)=1;
                end
            elseif n_x>32.5 && n_x<62.5 && n_y>=17 && n_y<=22
                dis2 = sqrt((n_x-c_x2)^2+(n_y-c_y2)^2);
                if dis2>=r
                    Macro_Label(xi,yi,:)=1;
                end       
            elseif n_x>117.5 && n_x<147.5 && n_y>=2 && n_y<=7
                dis3 = sqrt((n_x-c_x3)^2+(n_y-c_y3)^2);
                if dis3>=r
                    Macro_Label(xi,yi,:)=1;
                end
            elseif n_x>117.5 && n_x<147.5 && n_y>=17 && n_y<=22
                dis4 = sqrt((n_x-c_x4)^2+(n_y-c_y4)^2);
                if dis4>=r
                    Macro_Label(xi,yi,:)=1;
                end
            end
        end
    end

end

function Marco_Vtk(result, filepath, filename, varianame, thickness_ratio)
    [X, Y, Z] = size(result);
    num_element = X * Y * Z;
    file = fopen(strcat(filepath, filename, '.vtk'), 'w');
    fprintf(file, '%s\n', '# vtk DataFile Version 3.0');
    fprintf(file, '%s\n', 'This is used for visualization with Paraview');
    fprintf(file, '%s\n', 'ASCII');
    fprintf(file, '%s\n', 'DATASET RECTILINEAR_GRID');
    fprintf(file, 'DIMENSIONS %i %i %i\n', X + 1, Y + 1, Z + 1);
    fprintf(file, 'X_COORDINATES %i float\n', X + 1);
    X_coord = (0:X);
    fprintf(file, '%f ', X_coord);
    fprintf(file, '\n');
    fprintf(file, 'Y_COORDINATES %i float\n', Y + 1);
    Y_coord = (0:Y);
    fprintf(file, '%f ', Y_coord);
    fprintf(file, '\n');
    fprintf(file, 'Z_COORDINATES %i float\n', Z + 1);
    Z_coord = (0:Z)*thickness_ratio;
    fprintf(file, '%f ', Z_coord);
    fprintf(file, '\n');
    fprintf(file, '%s %i\n', 'CELL_DATA', num_element);
    fprintf(file, '%s %s %s\n', 'SCALARS', varianame, 'FLOAT');
    fprintf(file, '%s\n', 'LOOKUP_TABLE default');
    Phase = reshape(result, [], 1);
    fprintf(file, '%f\n', Phase);
    fclose(file);
end

% Calculate dog bones - simulate uniform Neumann boundary conditions

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

% Set the working directory
RVE_folder_load = strcat('FFT Dog');
if ~exist(RVE_folder_load,'dir')
    mkdir(RVE_folder_load); 
end
RVE_file_path_load = strcat(pwd,'\',RVE_folder_load,'\');

% 可视化初始模型信息
Marco_Vtk(RVE_Label,RVE_file_path_load,'Label','Label',RVE_L3/RVE_N3/(RVE_L1/RVE_N1));

% 赋予材料
Fiber_E22 = 300000.0; Fiber_Posi_rito_12 = 0.3;
Matrix_E = 0.001; Matrix_Posi_rito = 0.3;

Fiber_Lame = Fiber_E22*Fiber_Posi_rito_12/((1+Fiber_Posi_rito_12)*(1-2*Fiber_Posi_rito_12));
Fiber_G = Fiber_E22/(2.0*(1+Fiber_Posi_rito_12));
Matrix_Lame = Matrix_E*Matrix_Posi_rito/((1+Matrix_Posi_rito)*(1-2*Matrix_Posi_rito));
Matrix_G = Matrix_E/(2.0*(1+Matrix_Posi_rito));

RVE_Lame0 = (Fiber_Lame+Matrix_Lame)/2; RVE_G0 = (Fiber_G+Matrix_G)/2;

RVE_Lame_lib = Fiber_Lame*unit_matrix; RVE_Lame_lib(RVE_Label==0) = Matrix_Lame; RVE_Lame_lib(RVE_Label==2) = Matrix_Lame;
RVE_G_lib = Fiber_G*unit_matrix; RVE_G_lib(RVE_Label==0) = Matrix_G; RVE_G_lib(RVE_Label==2) = Matrix_G;


% 计算投影算子
RVE_ksi = 2*pi*[0:RVE_N1/2-1,-RVE_N1/2:-1]/RVE_L1;
RVE_eta = 2*pi*[0:RVE_N2/2-1,-RVE_N2/2:-1]/RVE_L2;
RVE_zdc = 2*pi*[0:RVE_N3/2-1,-RVE_N3/2:-1]/RVE_L3;

[RVE_G1111,RVE_G1122,RVE_G1133,RVE_G1112,RVE_G1113,RVE_G1123,...
  RVE_G2211,RVE_G2222,RVE_G2233,RVE_G2212,RVE_G2213,RVE_G2223,...
  RVE_G3311,RVE_G3322,RVE_G3333,RVE_G3312,RVE_G3313,RVE_G3323,...
  RVE_G1211,RVE_G1222,RVE_G1233,RVE_G1212,RVE_G1213,RVE_G1223,...
  RVE_G1311,RVE_G1322,RVE_G1333,RVE_G1312,RVE_G1313,RVE_G1323,...
  RVE_G2311,RVE_G2322,RVE_G2333,RVE_G2312,RVE_G2313,RVE_G2323] = ...
     Green_operateur_3D(RVE_L1,RVE_L2,RVE_L3,RVE_N1,RVE_N2,RVE_N3,...
     RVE_Lame0,RVE_G0,RVE_ksi,RVE_eta,RVE_zdc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Start FFT calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RVE_Load = [0,0,0,0,0,0]';

[RVE_sig11,RVE_sig22,RVE_sig33,RVE_sig23,RVE_sig13,RVE_sig12,...
    RVE_e11,RVE_e22,RVE_e33,RVE_e23,RVE_e13,RVE_e12] ...
    = FFT_Solver_3D(RVE_Load,RVE_N1,RVE_N2,RVE_N3,...
        RVE_Lame_lib,RVE_G_lib,RVE_Label,...
        RVE_G1111,RVE_G1122,RVE_G1133,RVE_G1112,RVE_G1113,RVE_G1123,...
        RVE_G2211,RVE_G2222,RVE_G2233,RVE_G2212,RVE_G2213,RVE_G2223,...
        RVE_G3311,RVE_G3322,RVE_G3333,RVE_G3312,RVE_G3313,RVE_G3323,...
        RVE_G1211,RVE_G1222,RVE_G1233,RVE_G1212,RVE_G1213,RVE_G1223,...
        RVE_G1311,RVE_G1322,RVE_G1333,RVE_G1312,RVE_G1313,RVE_G1323,...
        RVE_G2311,RVE_G2322,RVE_G2333,RVE_G2312,RVE_G2313,RVE_G2323);

figure(32);
aaa = RVE_sig11; aaa(RVE_Label~=1)=0;
imagesc(squeeze(aaa(:,:,1)));
colormap(jet);
colorbar;
set(gca, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none');

save("FFT_sig11.mat","RVE_sig11");
save("FFT_sig22.mat","RVE_sig22");
save("FFT_sig12.mat","RVE_sig12");

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

function [RVE_G1111,RVE_G1122,RVE_G1133,RVE_G1112,RVE_G1113,RVE_G1123,...
          RVE_G2211,RVE_G2222,RVE_G2233,RVE_G2212,RVE_G2213,RVE_G2223,...
          RVE_G3311,RVE_G3322,RVE_G3333,RVE_G3312,RVE_G3313,RVE_G3323,...
          RVE_G1211,RVE_G1222,RVE_G1233,RVE_G1212,RVE_G1213,RVE_G1223,...
          RVE_G1311,RVE_G1322,RVE_G1333,RVE_G1312,RVE_G1313,RVE_G1323,...
          RVE_G2311,RVE_G2322,RVE_G2333,RVE_G2312,RVE_G2313,RVE_G2323] = ...
             Green_operateur_3D(RVE_L1,RVE_L2,RVE_L3,RVE_N1,RVE_N2,RVE_N3,...
             RVE_Lame0,RVE_G0,RVE_ksi,RVE_eta,RVE_zdc)

    zero_matrix = zeros(RVE_N1,RVE_N2,RVE_N3);

    RVE_G1111 = zero_matrix; RVE_G1122 = zero_matrix; RVE_G1133 = zero_matrix;
    RVE_G1112 = zero_matrix; RVE_G1113 = zero_matrix; RVE_G1123 = zero_matrix;
                             RVE_G2222 = zero_matrix; RVE_G2233 = zero_matrix;
    RVE_G2212 = zero_matrix; RVE_G2213 = zero_matrix; RVE_G2223 = zero_matrix;
                                                      RVE_G3333 = zero_matrix;
    RVE_G3312 = zero_matrix; RVE_G3313 = zero_matrix; RVE_G3323 = zero_matrix;
    
    RVE_G1212 = zero_matrix; RVE_G1213 = zero_matrix; RVE_G1223 = zero_matrix;
    
                             RVE_G1313 = zero_matrix; RVE_G1323 = zero_matrix;
    
                                                      RVE_G2323 = zero_matrix;


    c1=1/4/RVE_G0; c2=(RVE_Lame0+RVE_G0)/RVE_G0/(RVE_Lame0+2*RVE_G0);

    for i=1:RVE_N1
        for j=1:RVE_N2
            for k=1:RVE_N3
                pr1 = RVE_ksi(i); pr2 = RVE_eta(j); pr3 = RVE_zdc(k);
                
                r1 = (RVE_N1/(0.5*RVE_L1))*sin(0.5*pr1*RVE_L1/RVE_N1)*cos(0.5*pr2*RVE_L2/RVE_N2)*cos(0.5*pr3*RVE_L3/RVE_N3);
                r2 = (RVE_N2/(0.5*RVE_L2))*cos(0.5*pr1*RVE_L1/RVE_N1)*sin(0.5*pr2*RVE_L2/RVE_N2)*cos(0.5*pr3*RVE_L3/RVE_N3);
                r3 = (RVE_N3/(0.5*RVE_L3))*cos(0.5*pr1*RVE_L1/RVE_N1)*cos(0.5*pr2*RVE_L2/RVE_N2)*sin(0.5*pr3*RVE_L3/RVE_N3);

                % r1 = pr1;
                % r2 = pr2;
                % r3 = pr3;

                d2=r1*r1+r2*r2+r3*r3;
                d4=d2*d2;

                if d2>0
                    RVE_G1111(i,j,k)=c1*4*r1*r1/d2-c2*(r1^4)/d4;
                    RVE_G1122(i,j,k)=-c2*(r1*r2)^2/d4;
                    RVE_G1133(i,j,k)=-c2*(r1*r3)^2/d4;
                    RVE_G1112(i,j,k)=c1*2*r1*r2/d2-c2*(r1^3)*r2/d4;
                    RVE_G1113(i,j,k)=c1*2*r1*r3/d2-c2*r3*(r1^3)/d4;
                    RVE_G1123(i,j,k)=-c2*r1*r1*r2*r3/d4;
                    
                    RVE_G2222(i,j,k)=c1*4*r2*r2/d2-c2*(r2^4)/d4;
                    RVE_G2233(i,j,k)=-c2*(r2*r3)^2/d4;
                    RVE_G2212(i,j,k)=c1*2*r1*r2/d2-c2*r1*(r2^3)/d4;
                    RVE_G2213(i,j,k)=-c2*r1*r2*r2*r3/d4;
                    RVE_G2223(i,j,k)=c1*2*r2*r3/d2-c2*(r2^3)*r3/d4;
                    
                    RVE_G3333(i,j,k)=c1*4*r3*r3/d2-c2*(r3^4)/d4; 
                    RVE_G3312(i,j,k)=-c2*r1*r2*r3*r3/d4;
                    RVE_G3313(i,j,k)=c1*2*r1*r3/d2-c2*r1*(r3^3)/d4;           
                    RVE_G3323(i,j,k)=c1*2*r2*r3/d2-c2*r2*r3^3/d4;   
                    
                    RVE_G1212(i,j,k)=c1*(r1*r1+r2*r2)/d2-c2*r1*r1*r2*r2/d4;
                    RVE_G1213(i,j,k)=c1*r3*r2/d2-c2*r1^2*r2*r3/d4; 
                    RVE_G1223(i,j,k)=c1*r1*r3/d2-c2*r1*r2^2*r3/d4;
                    
                    RVE_G1313(i,j,k)=c1*(r1*r1+r3*r3)/d2-c2*r1*r1*r3*r3/d4;
                    RVE_G1323(i,j,k)=c1*r1*r2/d2-c2*r1*r2*r3*r3/d4;  
                    
                    RVE_G2323(i,j,k)=c1*(r2*r2+r3*r3)/d2-c2*r2*r2*r3*r3/d4; 
                end
            end
        end
    end

    RVE_G2211 = RVE_G1122; RVE_G3311 = RVE_G1133; RVE_G3322 = RVE_G2233;

    RVE_G1112 = RVE_G1112*2.0; RVE_G1211 = RVE_G1112;
    RVE_G2212 = RVE_G2212*2.0; RVE_G1222 = RVE_G2212;
    RVE_G3312 = RVE_G3312*2.0; RVE_G1233 = RVE_G3312;
    RVE_G1113 = RVE_G1113*2.0; RVE_G1311 = RVE_G1113;
    RVE_G2213 = RVE_G2213*2.0; RVE_G1322 = RVE_G2213;
    RVE_G3313 = RVE_G3313*2.0; RVE_G1333 = RVE_G3313;
    RVE_G1123 = RVE_G1123*2.0; RVE_G2311 = RVE_G1123;
    RVE_G2223 = RVE_G2223*2.0; RVE_G2322 = RVE_G2223;
    RVE_G3323 = RVE_G3323*2.0; RVE_G2333 = RVE_G3323;

    RVE_G1213 = RVE_G1213*4.0; RVE_G1312 = RVE_G1213;
    RVE_G1223 = RVE_G1223*4.0; RVE_G2312 = RVE_G1223;
    RVE_G1323 = RVE_G1323*4.0; RVE_G2313 = RVE_G1323;

    RVE_G1212 = RVE_G1212*4.0; RVE_G1313 = RVE_G1313*4.0; RVE_G2323 = RVE_G2323*4.0;

    RVE_G1111(1,1,1) = 0; RVE_G1122(1,1,1) = 0; RVE_G1133(1,1,1) = 0;
    RVE_G1112(1,1,1) = 0; RVE_G1113(1,1,1) = 0; RVE_G1123(1,1,1) = 0;
    RVE_G2211(1,1,1) = 0; RVE_G2222(1,1,1) = 0; RVE_G2233(1,1,1) = 0;
    RVE_G2212(1,1,1) = 0; RVE_G2213(1,1,1) = 0; RVE_G2223(1,1,1) = 0;
    RVE_G3311(1,1,1) = 0; RVE_G3322(1,1,1) = 0; RVE_G3333(1,1,1) = 0;
    RVE_G3312(1,1,1) = 0; RVE_G3313(1,1,1) = 0; RVE_G3323(1,1,1) = 0;
    RVE_G1211(1,1,1) = 0; RVE_G1222(1,1,1) = 0; RVE_G1233(1,1,1) = 0;
    RVE_G1212(1,1,1) = 0; RVE_G1213(1,1,1) = 0; RVE_G1223(1,1,1) = 0;
    RVE_G1311(1,1,1) = 0; RVE_G1322(1,1,1) = 0; RVE_G1333(1,1,1) = 0;
    RVE_G1312(1,1,1) = 0; RVE_G1313(1,1,1) = 0; RVE_G1323(1,1,1) = 0;
    RVE_G2311(1,1,1) = 0; RVE_G2322(1,1,1) = 0; RVE_G2333(1,1,1) = 0;
    RVE_G2312(1,1,1) = 0; RVE_G2313(1,1,1) = 0; RVE_G2323(1,1,1) = 0;

end

function [RVE_sig11,RVE_sig22,RVE_sig33,RVE_sig23,RVE_sig13,RVE_sig12,...
            RVE_e11,RVE_e22,RVE_e33,RVE_e23,RVE_e13,RVE_e12] ...
            = FFT_Solver_3D(RVE_Load,RVE_N1,RVE_N2,RVE_N3,...
                RVE_Lame_lib,RVE_G_lib,RVE_Label,...
                RVE_G1111,RVE_G1122,RVE_G1133,RVE_G1112,RVE_G1113,RVE_G1123,...
                RVE_G2211,RVE_G2222,RVE_G2233,RVE_G2212,RVE_G2213,RVE_G2223,...
                RVE_G3311,RVE_G3322,RVE_G3333,RVE_G3312,RVE_G3313,RVE_G3323,...
                RVE_G1211,RVE_G1222,RVE_G1233,RVE_G1212,RVE_G1213,RVE_G1223,...
                RVE_G1311,RVE_G1322,RVE_G1333,RVE_G1312,RVE_G1313,RVE_G1323,...
                RVE_G2311,RVE_G2322,RVE_G2333,RVE_G2312,RVE_G2313,RVE_G2323)
    
    Err_tol = 1e-8;

    unit_matrix = ones(RVE_N1,RVE_N2,RVE_N3);
    % 首先初始化应变
    RVE_e11 = unit_matrix.*RVE_Load(1); RVE_e22 = unit_matrix.*RVE_Load(2); RVE_e33 = unit_matrix.*RVE_Load(3);
    RVE_e23 = unit_matrix.*RVE_Load(4); RVE_e13 = unit_matrix.*RVE_Load(5); RVE_e12 = unit_matrix.*RVE_Load(6);

    % 计算初始的应力
    RVE_sig11 = (RVE_Lame_lib+2.0*RVE_G_lib).*RVE_e11+RVE_Lame_lib.*RVE_e22+RVE_Lame_lib.*RVE_e33;
    RVE_sig22 = RVE_Lame_lib.*RVE_e11+(RVE_Lame_lib+2.0*RVE_G_lib).*RVE_e22+RVE_Lame_lib.*RVE_e33;
    RVE_sig33 = RVE_Lame_lib.*RVE_e11+RVE_Lame_lib.*RVE_e22+(RVE_Lame_lib+2.0*RVE_G_lib).*RVE_e33;
    RVE_sig23 = RVE_G_lib.*RVE_e23; RVE_sig13 = RVE_G_lib.*RVE_e13; RVE_sig12 = RVE_G_lib.*RVE_e12;

    RVE_sig11(RVE_Label==0) = 0; RVE_sig22(RVE_Label==0) = 0; RVE_sig12(RVE_Label==0) = 0;
    RVE_sig11(RVE_Label==0) = 0; RVE_sig22(RVE_Label==2) = 0; RVE_sig12(RVE_Label==2) = 0;
    RVE_sig11(1:50,11:110,:) = 100.0; RVE_sig11(RVE_N1-49:RVE_N1,11:110,:) = 100.0; 

    % 迭代求解
    Err_now = +inf; iter = 0;

    while Err_now>=Err_tol

        iter = iter+1;
        old_RVE_sig11 = RVE_sig11; old_RVE_sig22 = RVE_sig22; old_RVE_sig33 = RVE_sig33;
        old_RVE_sig23 = RVE_sig23; old_RVE_sig13 = RVE_sig13; old_RVE_sig12 = RVE_sig12;
        % 对应变和应力进行FFT
        F_RVE_e11 = fftn(RVE_e11); F_RVE_e22 = fftn(RVE_e22); F_RVE_e33 = fftn(RVE_e33);
        F_RVE_e23 = fftn(RVE_e23); F_RVE_e13 = fftn(RVE_e13); F_RVE_e12 = fftn(RVE_e12);
        F_RVE_sig11 = fftn(RVE_sig11); F_RVE_sig22 = fftn(RVE_sig22); F_RVE_sig33 = fftn(RVE_sig33);
        F_RVE_sig23 = fftn(RVE_sig23); F_RVE_sig13 = fftn(RVE_sig13); F_RVE_sig12 = fftn(RVE_sig12);
        
        % 计算L-S方程
        F_RVE_e11 = F_RVE_e11 - (RVE_G1111.*F_RVE_sig11 + RVE_G1122.*F_RVE_sig22 + RVE_G1133.*F_RVE_sig33 + RVE_G1123.*F_RVE_sig23 + RVE_G1113.*F_RVE_sig13 + RVE_G1112.*F_RVE_sig12); F_RVE_e11(1,1,1) = RVE_Load(1)*RVE_N1*RVE_N2*RVE_N3;
        F_RVE_e22 = F_RVE_e22 - (RVE_G2211.*F_RVE_sig11 + RVE_G2222.*F_RVE_sig22 + RVE_G2233.*F_RVE_sig33 + RVE_G2223.*F_RVE_sig23 + RVE_G2213.*F_RVE_sig13 + RVE_G2212.*F_RVE_sig12); F_RVE_e22(1,1,1) = RVE_Load(2)*RVE_N1*RVE_N2*RVE_N3;
        F_RVE_e33 = F_RVE_e33 - (RVE_G3311.*F_RVE_sig11 + RVE_G3322.*F_RVE_sig22 + RVE_G3333.*F_RVE_sig33 + RVE_G3323.*F_RVE_sig23 + RVE_G3313.*F_RVE_sig13 + RVE_G3312.*F_RVE_sig12); F_RVE_e33(1,1,1) = RVE_Load(3)*RVE_N1*RVE_N2*RVE_N3;
        F_RVE_e23 = F_RVE_e23 - (RVE_G2311.*F_RVE_sig11 + RVE_G2322.*F_RVE_sig22 + RVE_G2333.*F_RVE_sig33 + RVE_G2323.*F_RVE_sig23 + RVE_G2313.*F_RVE_sig13 + RVE_G2312.*F_RVE_sig12); F_RVE_e23(1,1,1) = RVE_Load(4)*RVE_N1*RVE_N2*RVE_N3;
        F_RVE_e13 = F_RVE_e13 - (RVE_G1311.*F_RVE_sig11 + RVE_G1322.*F_RVE_sig22 + RVE_G1333.*F_RVE_sig33 + RVE_G1323.*F_RVE_sig23 + RVE_G1313.*F_RVE_sig13 + RVE_G1312.*F_RVE_sig12); F_RVE_e13(1,1,1) = RVE_Load(5)*RVE_N1*RVE_N2*RVE_N3;
        F_RVE_e12 = F_RVE_e12 - (RVE_G1211.*F_RVE_sig11 + RVE_G1222.*F_RVE_sig22 + RVE_G1233.*F_RVE_sig33 + RVE_G1223.*F_RVE_sig23 + RVE_G1213.*F_RVE_sig13 + RVE_G1212.*F_RVE_sig12); F_RVE_e12(1,1,1) = RVE_Load(6)*RVE_N1*RVE_N2*RVE_N3;

        % FFT逆变换
        RVE_e11 = real(ifftn((F_RVE_e11))); RVE_e22 = real(ifftn((F_RVE_e22))); RVE_e33 = real(ifftn((F_RVE_e33)));
        RVE_e23 = real(ifftn((F_RVE_e23))); RVE_e13 = real(ifftn((F_RVE_e13))); RVE_e12 = real(ifftn((F_RVE_e12)));
        
        % 更新应力
        RVE_sig11 = (RVE_Lame_lib+2.0*RVE_G_lib).*RVE_e11+RVE_Lame_lib.*RVE_e22+RVE_Lame_lib.*RVE_e33;
        RVE_sig22 = RVE_Lame_lib.*RVE_e11+(RVE_Lame_lib+2.0*RVE_G_lib).*RVE_e22+RVE_Lame_lib.*RVE_e33;
        RVE_sig33 = RVE_Lame_lib.*RVE_e11+RVE_Lame_lib.*RVE_e22+(RVE_Lame_lib+2.0*RVE_G_lib).*RVE_e33;
        RVE_sig23 = RVE_G_lib.*RVE_e23; RVE_sig13 = RVE_G_lib.*RVE_e13; RVE_sig12 = RVE_G_lib.*RVE_e12;

        RVE_sig11(RVE_Label==0) = 0; RVE_sig22(RVE_Label==0) = 0; RVE_sig12(RVE_Label==0) = 0;
        RVE_sig11(RVE_Label==2) = 0; RVE_sig22(RVE_Label==2) = 0; RVE_sig12(RVE_Label==2) = 0;
        RVE_sig11(1:50,11:110,:) = 100.0; RVE_sig11(RVE_N1-49:RVE_N1,11:110,:) = 100.0; 

        % 计算误差
        Err_now = (norm(RVE_sig11-old_RVE_sig11,"fro")+norm(RVE_sig22-old_RVE_sig22,"fro")+norm(RVE_sig33-old_RVE_sig33,"fro")+...
                   norm(RVE_sig23-old_RVE_sig23,"fro")+norm(RVE_sig13-old_RVE_sig13,"fro")+norm(RVE_sig12-old_RVE_sig12,"fro"))/6/(RVE_N1)/(RVE_N2)/(RVE_N3);
        disp(['Iteration: ',num2str(iter),', Error: ',num2str(Err_now)]);
    end

end


% This program extracts and plots the defect potential

% NOTE: need to add treatment for [output_format] = 3 to feed EPW
%       with [eimp_mode] = 5 or 6

clear all;

%% material

% mat = 'Si';
% defect_name = 'Sb_Si';
% ifrelax = '/fixed';
% % ifrelax = '/relaxed';
% defect_type = 'S';
% % this should be consistent with DFT calculation, and also 
% %      the value used in scattering rate calculation, unit [aB]
% alat0_defect = 10.1802578337;  % Si, calculated using pz-hgh
% alat0_epw = 10.18101857;   % should be pz-n-nc

% mat = 'GaAs';
% defect_name = 'Zn_Ga';
% ifrelax = '/relaxed';
% defect_type = 'S';
% alat0_defect = 10.444418;    % GaAs
% alat0_epw = 10.4444371213;

mat = 'PbTe';
defect_name = 'Bi_Pb';
% defect_name = 'I_Te';
ifrelax = '/relaxed';
defect_type = 'S';
alat0_defect = 11.887552;
alat0_epw = 11.8873;

% mat = 'PbSe';
% % defect_name = 'Bi_Pb';
% defect_name = 'I_Se';
% ifrelax = '/relaxed';
% defect_type = 'S';
% alat0_defect = 11.2552919562;
% alat0_epw = 11.2553;

% mat = 'ZrNiSn';
% defect_name = 'Cu_Ni';
% % ifrelax = '/fixed';
% ifrelax = '/relaxed';
% defect_type = 'S';
% alat0_defect = 11.6913545169;   % ZrNiSn
% alat0_epw = 11.687027;

% mat = 'TiNiSn';
% defect_name = 'Nb_Ti';
% ifrelax = '/fixed';
% % ifrelax = '/relaxed';
% defect_type = 'S';
% alat0_defect = 11.3391373385;
% alat0_epw = 11.3351718627;

% mat = 'ZrCoSb';
% defect_name = 'Sc_Zr';
% % ifrelax = '/fixed';
% ifrelax = '/relaxed';
% defect_type = 'S';
% alat0_defect = 11.6011973852;   % ZrCoSb
% alat0_epw = 11.6053;

% mat = 'ZrCoBi';
% defect_name = 'Sc_Zr';
% ifrelax = '/relaxed';
% defect_type = 'S';
% alat0_defect = 11.7847937429;
% alat0_epw = 11.789239;

% mat = 'TiCoSb';
% defect_name = 'Y_Ti';
% % ifrelax = '/fixed';
% ifrelax = '/relaxed';
% defect_type = 'S';
% alat0_defect = 11.2594253624;
% alat0_epw = 11.2639779267;

% mat = 'NbFeSb';
% defect_name = 'Ti_Nb';
% % ifrelax = '/fixed';
% ifrelax = '/relaxed';
% defect_type = 'S';
% alat0_defect = 11.34039671;   % NbFeSb
% alat0_epw = 11.343926;

% mat = 'TaFeSb';
% defect_name = 'Ti_Ta';
% % ifrelax = '/fixed';
% ifrelax = '/relaxed';
% defect_type = 'S';
% alat0_defect = 11.3495280115;
% alat0_epw = 11.3529315843;

% mat = 'NbCoSn';
% defect_name = 'Zr_Nb';
% ifrelax = '';
% defect_type = 'S';
% alat0_defect = 11.32927192;   % NbCoSn
% alat0_epw = 11.328429;

% mat = 'Mg3Sb2';
% defect_name = 'V_Mg';
% defect_type = 'S';
% alat0 = 1;   % angstrom unit

charge = 1;
defect_pos = [0.00 0.00 0.00];  % relative to alat0 - primitive unit cell lattice vector
% defect_pos = [0.25 0.25 0.25];
% defect_pos = [0.50 0.50 0.50];

% for PbTe
% defect_pos = [-0.50 0.50 0.50];

% defect_pos = [0.111111111 0.222222222333333 0.210627926333333]; % crystal coordinate

% reference config
ref_config = 2;     % 1: separately in each defect folder, GaAs
                    % 2: collected in a main folder called R\, Si, ZrNiSn, etc.

% =========================================================================
% real space mesh for wave function (to be consistent with evc_check_r.dat
% nfft_wvfc = 40;  % Si, GaAs
% nfft_wvfc = 60;  % half Heuslers
nfft_wvfc = 48;  % PbTe
% =========================================================================

% defect                    
celltype = 'Conventional';
charge_R = 0;   % charge used for reference calculation
cellsize = 2;

model_charge = charge;

% lattice - conventional unit cell
a1 = [1,0,0];
a2 = [0,1,0];
a3 = [0,0,1];
alat = alat0_defect*cellsize;

% constant
aB = 0.529177249;  % Angstrom
Ryd = 13.605693;    % eV

% control
ioutput = 2;            % output = 0: check long range subtraction
                        % output on 1: original mesh; 2: interpolated mesh
                        % output = 3: point-by-point check Coulomb potential subtraction
                        % output = -1: point-by-point, compare DFT and Coulomb potential
                        
%% ========================================================================
% setting
read_dft = 1;
if ioutput == -1
    point_by_point = 1;
    ivcoul = 0;
    icheck = 1;
elseif ioutput == 0
    point_by_point = 0;     % if 1 then do point_by_point correction of electrostatic potential
    ivcoul = 1;             % if 1 subtracts Coulomb model potential from total potential
    icheck = 1;
elseif ioutput == 1 || ioutput == 2
    point_by_point = 1;
    ivcoul = 0;
    icheck = 0;
else
    point_by_point = 1;
    ivcoul = 1;
    icheck = 1;    
end

%%
defect_type_R = 'R';
nfft_interpl = cellsize*nfft_wvfc;     % interpolated mesh for supercell, should match with nfft_wvfc

% new
if (ref_config == 1)
    file_ionhar_R = strcat(mat,'/',defect_name,'/q',num2str(charge),'/n',num2str(cellsize),'/',mat,'.',defect_type_R,'_q',num2str(charge_R),'_',celltype(1),num2str(cellsize),'.ionharpot.xsf');
    file_ion_R = strcat(mat,'/',defect_name,'/q',num2str(charge),'/n',num2str(cellsize),'/',mat,'.',defect_type_R,'_q',num2str(charge_R),'_',celltype(1),num2str(cellsize),'.ionpot.xsf');
    file_tot_R = strcat(mat,'/',defect_name,'/q',num2str(charge),'/n',num2str(cellsize),'/',mat,'.',defect_type_R,'_q',num2str(charge_R),'_',celltype(1),num2str(cellsize),'.totpot.xsf');
else
    file_ionhar_R = strcat(mat,'/R',ifrelax,'/q',num2str(charge_R),'/n',num2str(cellsize),'/',mat,'.',defect_type_R,'_q',num2str(charge_R),'_',celltype(1),num2str(cellsize),'.ionharpot.xsf');
    file_ion_R = strcat(mat,'/R',ifrelax,'/q',num2str(charge_R),'/n',num2str(cellsize),'/',mat,'.',defect_type_R,'_q',num2str(charge_R),'_',celltype(1),num2str(cellsize),'.ionpot.xsf');
    file_tot_R = strcat(mat,'/R',ifrelax,'/q',num2str(charge_R),'/n',num2str(cellsize),'/',mat,'.',defect_type_R,'_q',num2str(charge_R),'_',celltype(1),num2str(cellsize),'.totpot.xsf');    
end
file_ionhar_S = strcat(mat,'/',defect_name,ifrelax,'/q',num2str(charge),'/n',num2str(cellsize),'/',mat,'.',defect_type,'_',defect_name,'_q',num2str(charge),'_',celltype(1),num2str(cellsize),'.ionharpot.xsf');
file_ion_S = strcat(mat,'/',defect_name,ifrelax,'/q',num2str(charge),'/n',num2str(cellsize),'/',mat,'.',defect_type,'_',defect_name,'_q',num2str(charge),'_',celltype(1),num2str(cellsize),'.ionpot.xsf');
file_tot_S = strcat(mat,'/',defect_name,ifrelax,'/q',num2str(charge),'/n',num2str(cellsize),'/',mat,'.',defect_type,'_',defect_name,'_q',num2str(charge),'_',celltype(1),num2str(cellsize),'.totpot.xsf');

if (read_dft == 1)
    %% extract ionic + Hartree potential
    size6 = [6 inf];

    iu_R = fopen(file_ionhar_R,'r');
    iu_S = fopen(file_ionhar_S,'r');
    s = fgetl(iu_R);
    while (length(strfind(s,'DATAGRID_3D_UNKNOWN')) == 0)
        s = fgetl(iu_R);     
    end
    s = fgetl(iu_S);
    while (length(strfind(s,'DATAGRID_3D_UNKNOWN')) == 0)
        s = fgetl(iu_S);
    end

    s = fgetl(iu_R);
    scan_l = sscanf(s,'%d %d %d');
    nfftx = scan_l(1);
    nffty = scan_l(2);
    nfftz = scan_l(3);
    nfft = nfftx*nffty*nfftz;

    v_ionhar_R = zeros(nfft,1);
    v_ionhar_S = zeros(nfft,1);

    s = fgetl(iu_R);
    at = zeros(3,3);

    s = fgetl(iu_R);
    scan_l = sscanf(s,'%f %f %f');
    at(1:3,1) = scan_l(1:3);    % at in [Angstrom]
    s = fgetl(iu_R);
    scan_l = sscanf(s,'%f %f %f');
    at(1:3,2) = scan_l(1:3);
    s = fgetl(iu_R);
    scan_l = sscanf(s,'%f %f %f');
    at(1:3,3) = scan_l(1:3);

    s = fgetl(iu_S);
    s = fgetl(iu_S);
    s = fgetl(iu_S);
    s = fgetl(iu_S);
    s = fgetl(iu_S);

    pot_ = fscanf(iu_R,'%f %f %f %f %f %f',size6);
    v_ionhar_R = pot_';
    pot_ = fscanf(iu_S,'%f %f %f %f %f %f',size6);
    v_ionhar_S = pot_';

    fclose(iu_R);
    fclose(iu_S);

    % --- perturbed potential
    dv_ionhar_ = v_ionhar_S - v_ionhar_R;
    nrow = length(dv_ionhar_);
    dv_ionhar = zeros(1,nfft);
    for i = 1:6;
        if (6*(nrow-1)+i <= nfft)
            dv_ionhar(i:6:6*(nrow-1)+i) = dv_ionhar_(1:nrow,i);
        else
            dv_ionhar(i:6:6*(nrow-1)) = dv_ionhar_(1:nrow-1,i);
        end
    end

    %% extract ionic potential
    iu_R = fopen(file_ion_R,'r');
    iu_S = fopen(file_ion_S,'r');
    s = fgetl(iu_R);           
    while (length(strfind(s,'DATAGRID_3D_UNKNOWN')) == 0)
        s = fgetl(iu_R);     
    end
    s = fgetl(iu_S);
    while (length(strfind(s,'DATAGRID_3D_UNKNOWN')) == 0)
        s = fgetl(iu_S);
    end

    s = fgetl(iu_R);
    s = fgetl(iu_R);
    s = fgetl(iu_R);
    s = fgetl(iu_R);
    s = fgetl(iu_R);
    s = fgetl(iu_S);
    s = fgetl(iu_S);
    s = fgetl(iu_S);
    s = fgetl(iu_S);
    s = fgetl(iu_S);

    v_ion_R = zeros(nfft,1);
    v_ion_S = zeros(nfft,1);
    pot_ = fscanf(iu_R,'%f %f %f %f %f %f',size6);
    v_ion_R = pot_';
    pot_ = fscanf(iu_S,'%f %f %f %f %f %f',size6);
    v_ion_S = pot_';

    fclose(iu_R);
    fclose(iu_S);

    % --- perturbed potential
    dv_ion_ = v_ion_S - v_ion_R;
    nrow = length(dv_ion_);
    dv_ion = zeros(1,nfft);
    for i = 1:6;
        if (6*(nrow-1)+i <= nfft)
            dv_ion(i:6:6*(nrow-1)+i) = dv_ion_(1:nrow,i);
        else
            dv_ion(i:6:6*(nrow-1)) = dv_ion_(1:nrow-1,i);
        end
    end

    %% extract total (ionic + Hartree + xc) potential
    iu_R = fopen(file_tot_R,'r');
    iu_S = fopen(file_tot_S,'r');
    s = fgetl(iu_R);           
    while (length(strfind(s,'DATAGRID_3D_UNKNOWN')) == 0)
        s = fgetl(iu_R);     
    end
    s = fgetl(iu_S);
    while (length(strfind(s,'DATAGRID_3D_UNKNOWN')) == 0)
        s = fgetl(iu_S);
    end

    s = fgetl(iu_R);
    s = fgetl(iu_R);
    s = fgetl(iu_R);
    s = fgetl(iu_R);
    s = fgetl(iu_R);
    s = fgetl(iu_S);
    s = fgetl(iu_S);
    s = fgetl(iu_S);
    s = fgetl(iu_S);
    s = fgetl(iu_S);

    v_tot_R = zeros(nfft,1);
    v_tot_S = zeros(nfft,1);
    pot_ = fscanf(iu_R,'%f %f %f %f %f %f',size6);
    v_tot_R = pot_';
    pot_ = fscanf(iu_S,'%f %f %f %f %f %f',size6);
    v_tot_S = pot_';

    fclose(iu_R);
    fclose(iu_S);

    % --- perturbed potential
    dv_tot_ = v_tot_S - v_tot_R;
    nrow = length(dv_tot_);
    dv_tot = zeros(1,nfft);
    for i = 1:6
        if (6*(nrow-1)+i <= nfft)
            dv_tot(i:6:6*(nrow-1)+i) = dv_tot_(1:nrow,i);
        else
            dv_tot(i:6:6*(nrow-1)) = dv_tot_(1:nrow-1,i);
        end
    end

    %  ---
    dv_har = dv_ionhar - dv_ion;
    dv_xc = dv_tot - dv_ionhar;

    % unit correction, to [eV]
    dv_tot = dv_tot * Ryd;
    dv_ion = dv_ion * Ryd;
    dv_har = dv_har * Ryd;
    dv_xc  = dv_xc  * Ryd;
end


%% interpolation to commensurate mesh
% -- original mesh
i0 = [1:nfft];
j0 = [1:nfft];
k0 = [1:nfft];
i0 = mod(i0-1,nfftx);
j0 = mod(floor((j0-1)/(nfftx)),nffty);
k0 = mod(floor((k0-1)/(nfftx)/(nffty)),nfftz);

x0 = (i0/(nfftx-1))*at(1,1) + (j0/(nffty-1))*at(1,2) + (k0/(nfftz-1))*at(1,3);
y0 = (i0/(nfftx-1))*at(2,1) + (j0/(nffty-1))*at(2,2) + (k0/(nfftz-1))*at(2,3);
z0 = (i0/(nfftx-1))*at(3,1) + (j0/(nffty-1))*at(3,2) + (k0/(nfftz-1))*at(3,3);

dv_tot_3d = zeros(nfftx,nffty,nfftz);
for i = 1:nfft
    dv_tot_3d(i0(i)+1,j0(i)+1,k0(i)+1) = dv_tot(i);
end

% -- interpolated mesh
i0_int = [1:nfft_interpl^3];
j0_int = [1:nfft_interpl^3];
k0_int = [1:nfft_interpl^3];
i0_int = mod(i0_int-1,nfft_interpl);
j0_int = mod(floor((j0_int-1)/(nfft_interpl)),nfft_interpl);
k0_int = mod(floor((k0_int-1)/(nfft_interpl)/(nfft_interpl)),nfft_interpl);

x0_int = (i0_int/nfft_interpl)*at(1,1) + (j0_int/nfft_interpl)*at(1,2) + (k0_int/nfft_interpl)*at(1,3);
y0_int = (i0_int/nfft_interpl)*at(2,1) + (j0_int/nfft_interpl)*at(2,2) + (k0_int/nfft_interpl)*at(2,3);
z0_int = (i0_int/nfft_interpl)*at(3,1) + (j0_int/nfft_interpl)*at(3,2) + (k0_int/nfft_interpl)*at(3,3);

% interpolating cubic corner
i1 = floor(i0_int/nfft_interpl*(nfftx-1)) + 1;
j1 = floor(j0_int/nfft_interpl*(nffty-1)) + 1;
k1 = floor(k0_int/nfft_interpl*(nfftz-1)) + 1;
ratio_i = (i0_int/nfft_interpl*(nfftx-1)) - (i1-1);
ratio_j = (j0_int/nfft_interpl*(nffty-1)) - (j1-1);
ratio_k = (k0_int/nfft_interpl*(nfftz-1)) - (k1-1);

for i = 1:nfft_interpl^3
    dv_tot_interpl(i) = dv_tot_3d(i1(i),  j1(i),  k1(i))   * (1-ratio_i(i))*(1-ratio_j(i))*(1-ratio_k(i)) + ...
                        dv_tot_3d(i1(i),  j1(i),  k1(i)+1) * (1-ratio_i(i))*(1-ratio_j(i))*   ratio_k(i)  + ...
                        dv_tot_3d(i1(i),  j1(i)+1,k1(i))   * (1-ratio_i(i))*   ratio_j(i) *(1-ratio_k(i)) + ...
                        dv_tot_3d(i1(i),  j1(i)+1,k1(i)+1) * (1-ratio_i(i))*   ratio_j(i) *   ratio_k(i)  + ...
                        dv_tot_3d(i1(i)+1,j1(i),  k1(i))   *    ratio_i(i) *(1-ratio_j(i))*(1-ratio_k(i)) + ...
                        dv_tot_3d(i1(i)+1,j1(i),  k1(i)+1) *    ratio_i(i) *(1-ratio_j(i))*   ratio_k(i)  + ...
                        dv_tot_3d(i1(i)+1,j1(i)+1,k1(i))   *    ratio_i(i) *   ratio_j(i) *(1-ratio_k(i)) + ...
                        dv_tot_3d(i1(i)+1,j1(i)+1,k1(i)+1) *    ratio_i(i) *   ratio_j(i) *   ratio_k(i);
end



%% find distance to closet impurity and correct

rshift = zeros(3,8);
icopy = 0;
for i = 1:2
    for j = 1:2
        for k = 1:2
            icopy = icopy + 1;
            rshift(:,icopy) = (i-1)*at(:,1) + (j-1)*at(:,2) + (k-1)*at(:,3);
        end
    end
end

xd = defect_pos(1) * alat0_defect*aB;
yd = defect_pos(2) * alat0_defect*aB;
zd = defect_pos(3) * alat0_defect*aB;

% -- potential on original mesh
nlist = 0;
rdist = zeros(1,8);
rp = zeros(1,1);
for i = 1:nfft
    i0 = mod(i-1,nfftx);
    j0 = mod(floor((i-1)/(nfftx)),nffty);
    k0 = mod(floor((i-1)/(nfftx)/(nffty)),nfftz);    

    if ((i0 == nfftx-1) || (j0 == nffty-1) || (k0 == nfftz-1)) 
        % remove repeated points
        continue
    end
    for j = 1:8
        rdist(j) = ((x0(i)-xd-rshift(1,j)).^2 + (y0(i)-yd-rshift(2,j)).^2 + (z0(i)-zd-rshift(3,j)).^2).^0.5;
    end
    rdist_min_ = min(rdist);
    % =============================================
    % the rdist_min found below is not uniquely set, should check if
    % duplicates have any effect
    % =============================================
    for j = 1:8
        % make sure the area is symmetric with respect to the origin
        if (rdist(j) == rdist_min_)
            nlist = nlist + 1;
            rp(nlist) = rdist(j);
            xc(nlist) = x0(i) - xd - rshift(1,j);
            yc(nlist) = y0(i) - yd - rshift(2,j);
            zc(nlist) = z0(i) - zd - rshift(3,j);
            dv_tot_c(nlist) = dv_tot(i);
        end
    end
end

% -- potential on interpolated mesh
nlist_int = 0;
rdist = zeros(1,8);
rp_int = zeros(1,1);
for i = 1:nfft_interpl^3
    i0 = mod(i-1,nfft_interpl);
    j0 = mod(floor((i-1)/(nfft_interpl)),nfft_interpl);
    k0 = mod(floor((i-1)/(nfft_interpl)/(nfft_interpl)),nfft_interpl);

    for j = 1:8
        rdist(j) = ((x0_int(i)-xd-rshift(1,j)).^2 + (y0_int(i)-yd-rshift(2,j)).^2 + (z0_int(i)-zd-rshift(3,j)).^2).^0.5;
    end
    rdist_min_ = min(rdist);
    for j = 1:8
        % make sure the area is symmetric with respect to the origin
        if (rdist(j) == rdist_min_)
            nlist_int = nlist_int + 1;
            rp_int(nlist_int) = rdist(j);
            xc_int(nlist_int) = x0_int(i) - xd - rshift(1,j);
            yc_int(nlist_int) = y0_int(i) - yd - rshift(2,j);
            zc_int(nlist_int) = z0_int(i) - zd - rshift(3,j);
            dv_tot_int_c(nlist_int) = dv_tot_interpl(i);
            break;
        end
    end
end

% unit conversion - in unit of primitive unit cell lattice constant
xc = xc / alat0_defect / aB;
yc = yc / alat0_defect / aB;
zc = zc / alat0_defect / aB;
xc_int = xc_int / alat0_defect / aB;
yc_int = yc_int / alat0_defect / aB;
zc_int = zc_int / alat0_defect / aB;

% convert to [Ryd]
dv_tot_c = dv_tot_c / Ryd;
dv_tot_int_c = dv_tot_int_c / Ryd;


%% output DFT defect potential to file

if (ioutput == 1)
    %
    output = fopen('dv_tot.dat','w');
    fprintf(output,'%12.8f\n',model_charge);
    fprintf(output,' %d %d %d\n',nfftx-1,nffty-1,nfftz-1);
    fprintf(output,'%12.8f\n',alat0_epw);
    fprintf(output,'%12.8f %12.8f %12.8f\n',defect_pos(:));    
    fprintf(output,'%12.8f %12.8f %12.8f\n',at(:,1)/alat0_defect/aB);
    fprintf(output,'%12.8f %12.8f %12.8f\n',at(:,2)/alat0_defect/aB);
    fprintf(output,'%12.8f %12.8f %12.8f\n',at(:,3)/alat0_defect/aB);
    fprintf(output,' %d\n',nlist);
    for i = 1:nlist
        i
        fprintf(output,'%12.8f %12.8f %12.8f %12.8f\n',xc(i),yc(i),zc(i),dv_tot_c(i));
    end
    fclose(output);    
    %
elseif (ioutput == 2)
    %
    output = fopen('dv_tot.dat','w');
    fprintf(output,'%12.8f\n',model_charge);
    fprintf(output,' %d %d %d\n',nfft_interpl,nfft_interpl,nfft_interpl);
    fprintf(output,'%12.8f\n',alat0_epw);
    fprintf(output,'%12.8f %12.8f %12.8f\n',defect_pos(:));        
    fprintf(output,'%12.8f %12.8f %12.8f\n',at(:,1)/alat0_defect/aB);
    fprintf(output,'%12.8f %12.8f %12.8f\n',at(:,2)/alat0_defect/aB);
    fprintf(output,'%12.8f %12.8f %12.8f\n',at(:,3)/alat0_defect/aB);
    fprintf(output,' %d\n',nlist_int);
    for i = 1:nlist_int
        i
        fprintf(output,'%12.8f %12.8f %12.8f %12.8f\n',xc_int(i),yc_int(i),zc_int(i),dv_tot_int_c(i));
    end
    fclose(output);
    %
end


%% check long range potential
if (point_by_point == 1 && icheck == 1)
%  point-by-point electrostatic impurity potential    
    at1_ = at(:,1)' / alat / aB;
    at2_ = at(:,2)' / alat / aB;
    at3_ = at(:,3)' / alat / aB;

    index = 0;
    rmodel = zeros(1,nlist);
    vcoul = zeros(1,nlist);
    if (ivcoul == 1)
        for i3 = 1:nfftz
%         for i3 = 50:60
            i3
            tic
            for i2 = 1:nffty
%             for i2 = 50:60
                for i1 = 1:nfftx
%                 for i1 = 50:60
%                     nlist = i1 + (i2-1)*nfftx + (i3-1)*nfftx*nffty;
                    if ((i1 == nfftx) || (i2 == nffty) || (i3 == nfftz)) 
                        % jump over repeated points
                        continue
                    end
                    index = index + 1;
                    rp_ = [xc(index) yc(index) zc(index)]'*(alat0_defect/alat);  % need to normalize to the supercell lattice
%                     rp = ((i1-1)/(nfftx-1))*at1_ + ((i2-1)/(nffty-1))*at2_ + ((i3-1)/(nfftz-1))*at3_;
                    vcoul(index) = V_Coulomb(mat, model_charge, at1_, at2_, at3_, alat, rp_);
                end
            end
            toc
        end
    end

    dv_tot_sr = dv_tot_c*Ryd - vcoul;   % in energy unit [eV]

    % find potential far away
%     rdistmax = 0;
%     for i = 1:nfft
%         if (rdistmax < rdist_min(i))
%             i_rmax = i;
%             rdistmax = rdist_min(i);
%             dvtot_far = dv_tot(i);
%             vcoul_far = vcoul(i);
%             dvtotsr_far = dv_tot_sr(i);
%         end
%     end
    
% %     plannar-averaged potential based on point-by-point model potential
%     rmodel_pa = zeros(nn-1,1);
%     vcoul_pa = zeros(nn-1,1);
%     for i1 = 1:nn-1
%         rmodel_pa(i1) = (i1/nn)*norm(a1,2)*alat*aB;
%         for i2 = 1:nn
%             for i3 = 1:nn
%                 rp = (i1/nn)*a1 + (i2/nn)*a2 + (i3/nn)*a3;
%                 vcoul_pa(i1) = vcoul_pa(i1) + V_Coulomb(mat, model_charge, a1, a2, a3, alat, rp);
%             end
%         end
%         vcoul_pa(i1) = vcoul_pa(i1)/nn/nn;
%     end

elseif (point_by_point == 0 && icheck == 1)
%  planar averaged potential y-z plane

    rdist_pa = zeros(nfftx,1);
    for ia = 1:nfftx
        rdist_pa(ia) = ((ia-1)/(nfftx-1))*norm(at(:,1),2);
    end

    dv_tot_pa = zeros(nfftx,1);
    dv_ion_pa = zeros(nfftx,1);
    dv_har_pa = zeros(nfftx,1);
    dv_xc_pa = zeros(nfftx,1);
    for i = 1:nfft
        ia = mod(i-1,nfftx) + 1;
        dv_tot_pa(ia) = dv_tot_pa(ia) + dv_tot(i);
        dv_ion_pa(ia) = dv_ion_pa(ia) + dv_ion(i);
        dv_har_pa(ia) = dv_har_pa(ia) + dv_har(i);
        dv_xc_pa(ia)  = dv_xc_pa(ia)  + dv_xc(i);
    end
    dv_tot_pa = dv_tot_pa / nffty / nfftz;
    dv_ion_pa = dv_ion_pa / nffty / nfftz;
    dv_har_pa = dv_har_pa / nffty / nfftz;
    dv_xc_pa  = dv_xc_pa  / nffty / nfftz;

    % compute model potential along DFT points
    dv_pa_1D = zeros(nfftx,1);
    dv_pa_1D_sr = zeros(nfftx,1);
    if (ivcoul == 1)
        x_half = 0.5*norm(at(:,1),2)/alat/aB;
        x_avg = x_half + defect_pos(1)*alat0_defect/alat;
        for ia = 1:nfftx
            xp = rdist_pa(ia)/alat/aB - x_avg;  % in [alat]
            if (xp < -x_half)
                xp = xp + 2*x_half;
            end
            dv_pa_1D(ia) = V_Coulomb_1D(mat, model_charge, a1, a2, a3, alat, xp);
        end
    end
    
    dv_pa_1D_sr = dv_tot_pa - dv_pa_1D;
    % dv_pa_1D_sr = dv_har_pa + dv_ion_pa - dv_pa_1D;
end


%% plot

if (point_by_point == 1)
    figure(1);
    hold on;

%     plot(rp,dv_tot_c,'.','markersize',8);
%     plot(rp_int,dv_tot_int_c,'.','markersize',8);
    
    plot(rp,dv_tot_c*Ryd,'.','markersize',8);
%     plot(rp,dv_tot_sr,'.','markersize',8);
%     plot(rp,vcoul,'.','markersize',8);
    
    % compare DFT and Coulomb potential
%     plot(rp,dv_tot_c*Ryd,'.','markersize',8);
    
    % TiCoSb
%     dielec = 19.09;
    % TiNiSn
%     dielec = 22.92;
    % ZrNiSn
%     dielec = 20.90;
    % ZrCoSb
%     dielec = 17.87;
    % Si
    dielec = 11.68;
    
    rp_l = [0.01:0.01:10]/(1e10);
    V_coulomb = ((model_charge*(-1.602e-19))/(dielec*(8.854187e-12)))/4/pi./rp_l;
    plot(rp_l*(1e10),V_coulomb,'--','linewidth',2);
    
    % save the data
    dv_mat(:,1) = rp';
    dv_mat(:,2) = dv_tot_c'*Ryd;
    dv_coul(:,1) = rp_l'*(1e10);
    dv_coul(:,2) = V_coulomb';

    axis([0 10 -15 10]);
    box on;
    set(gca,'fontsize',16);
    xlabel('Distance (Angstrom)');
    ylabel('Potential energy (eV)');
    legend('Central cell potential','Coulomb potential','box','off');
    
    % % Ionic
    % figure(4);
    % 
    % plot(rdist_min,dv_tot,'.','markersize',5);
    % plot(rdist_min,dv_ion,'.','markersize',5);
    % 
    % % Hartreen
    % figure(5);
    % 
    % plot(rdist_min,dv_har,'.','markersize',5);
    % 
    % % exchange-correlation
    % figure(6);
    % 
    % plot(rdist_min,dv_xc,'.','markersize',5);

elseif (point_by_point == 0 && icheck == 1)
    figure;
    hold on;

    plot(rdist_pa,dv_tot_pa,'o-','markersize',5);
    % plot(rdist_pa,dv_ion_pa,'o-','markersize',5);
    % plot(rdist_pa,dv_har_pa,'o-','markersize',5);
    % plot(rdist_pa,dv_har_pa+dv_ion_pa,'o-','markersize',5);
    % plot(rdist_pa,dv_xc_pa,'o-','markersize',5);

    plot(rdist_pa,dv_pa_1D,'o-.','markersize',5);
    plot(rdist_pa,dv_pa_1D_sr,'o--','markersize',5);
    
end


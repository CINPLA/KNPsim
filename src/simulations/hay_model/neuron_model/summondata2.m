function summondata2
% Returns
% N = number of voxels
% times = times
% jk,   jna, jca: Ion fluxes onto voxels (N x times)
% jx:   flux of unspecified ion into voxels (assumed valence -1)
%       jx contains il, ih and isyn
% icap: capacitive current into voxel
% isyn: synaptic current into voxel (contained in jx)

N = 196
iica = [];
iina = [];
iik = [];
iix = [];
iicap = [];
it = [];
iimemb = [];
iisyn = [];

str1 = 'new_currsums_parts_10000areagsynsmediumtau_fixeddt_type2_amp4.2e-05_tstop10000.0_nseg20_dt0.025_';
str2 = '_comb_summed.mat';

for i = 1:1
strseed = ['seed',num2str(i)];
strfile = [str1,strseed,str2];
load(strfile);

iion0 = ik + ina + ica + il + ih;
isyn = imemb - (iion0 + icap);
ix = il + ih + isyn; % all currents of unspecified ion species x

ts2 = ts(:,ts>=1600)-1600 + 8400*(i-1);
ica = ica(:,ts>=1600);
ina = ina(:,ts>=1600);
ik = ik(:,ts>=1600);
ix = ix(:,ts>=1600);
icap = icap(:,ts>=1600);
imemb = imemb(:,ts>=1600);
isyn = isyn(:,ts>=1600);

it = [it,ts2];
iica = [iica,ica];
iina = [iina, ina];
iik = [iik,ik];
iix = [iix,ix];
iicap = [iicap,icap];
iimemb = [iimemb,imemb];
iisyn = [iisyn,isyn];

end

% Convert to SI units
ts = it/1000; %s
icap = iicap*1e-9; % Ampere
imemb = iimemb*1e-9;
isyn = iisyn*1e-9;
ik = iik*1e-9; ina = iina*1e-9; ica = iica*1e-9; ix = iix*1e-9;

% Convert currents to fluxes (j_k = i_k/(z_k*F))
F = 96485.3365; % C/mol
jk = ik/F; % mol/s
jna = ina/F;
jca = ica/2/F;
jx = -ix/F; % I define x as unknown negative ion with valence 1
x = x*1e-6; % meters
y = y*1e-6; % meters
z = z*1e-6; % meters


save('revdata2.mat','N', 'ts', 'jk', 'jna', 'jca', 'jx', 'icap', 'isyn', 'imemb', 'x', 'y', 'z');
end


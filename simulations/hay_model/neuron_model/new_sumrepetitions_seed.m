% combinemattomat_fixeddt.m
% A MATLAB script for summing the  the files saved by running "calcsumcurr_manyareagsynmediumtau_parts_fixeddt.py 20 0.025 0.000042 10000 10000 2 myseed 200"
% Expects that the variable "iseed" has been initialized. iseed=1 corresponds to myseed=1,...,10 in combinemattomat_fixeddt.m, iseed=2 corresponds to
% myseed=11,...,20, etc.
% Tuomo Maki-Marttunen, 2014-2016

iseed = 1

synloctype = 2;

dt = 0.1;
T=10000;
Nts = round(T/dt)-1;
ts = dt*(0:Nts-1);
Nvox=196;


ica = zeros(Nvox,Nts);
icap = zeros(Nvox,Nts);
il = zeros(Nvox,Nts);
ik = zeros(Nvox,Nts);
ih = zeros(Nvox,Nts);
ina = zeros(Nvox,Nts);
imemb = zeros(Nvox,Nts);
VtimesA = zeros(Nvox,Nts);

Vsoma = cell(1,10);
ts_syn = cell(1,10);
part_syn = cell(1,10);
for irep=1:1
  disp(['Loading ' num2str(iseed)]);
  A=load(['new_currsums_parts_10000areagsynsmediumtau_fixeddt_type' num2str(synloctype) '_amp4.2e-05_tstop10000.0_nseg20_dt0.025_seed' num2str((iseed-1)*10+irep) '_comb200.0.mat']);
  ica = ica + interpolate_multidim(A.times,A.ica,ts);  
  icap = icap + interpolate_multidim(A.times,A.icap,ts);  
  il = il + interpolate_multidim(A.times,A.il,ts);
  ik = ik + interpolate_multidim(A.times,A.ik,ts);
  ih = ih + interpolate_multidim(A.times,A.ih,ts);
  ina = ina + interpolate_multidim(A.times,A.ina,ts);
  imemb = imemb + interpolate_multidim(A.times,A.imemb,ts);
  VtimesA = VtimesA + interpolate_multidim(A.times,A.VtimesA,ts);
  Vsoma{irep} = interpolate(A.times,A.Vsoma,ts);
  ts_syn{irep} = A.ts_syn;
  part_syn{irep} = A.part_syn;
  x = A.x; 
  y = A.y; 
  z = A.z;
end
save(['new_currsums_parts_10000areagsynsmediumtau_fixeddt_type' num2str(synloctype) '_amp4.2e-05_tstop10000.0_nseg20_dt0.025_seed' num2str(iseed) '_comb_summed.mat'],'ts','ica','icap','il','ik','ih','ina','imemb','VtimesA','Vsoma','ts_syn','part_syn', 'x', 'y', 'z');


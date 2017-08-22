%function to save F, P, SO, SN with a given file number
function  tmp = Save_FPIONAGfp(Fn, Pn, In, SOn, SNn, SAn, SGn, TBn, TTn, theta_f_n, theta_p_n, imax, file_NO)

global NX NY Ind_B Ind_T

F2D = reshape(Fn,NX+1,NY+1);
F2D = F2D';

P2D = reshape(Pn,NX+1,NY+1);
P2D = P2D';

I2D = reshape(In,NX+1,NY+1);
I2D = I2D';

SO2D = reshape(SOn,NX+1,NY+1);
SO2D = SO2D';

SN2D = reshape(SNn,NX+1,NY+1);
SN2D = SN2D';

SA2D = reshape(SAn,NX+1,NY+1);
SA2D = SA2D';

SG2D = reshape(SGn,NX+1,NY+1);
SG2D = SG2D';

if Ind_B == 1
    TB2D = reshape(TBn,NX+1,NY+1);
    TB2D = TB2D';
end

if Ind_T == 1
    TT2D = reshape(TTn,NX+1,NY+1);
    TT2D = TT2D';
end

theta_f_2D = reshape(theta_f_n,NX+1,NY+1);
theta_f_2D = theta_f_2D';

theta_p_2D = reshape(theta_p_n,NX+1,NY+1);
theta_p_2D = theta_p_2D';

fname_F = sprintf('%s_%d','./data/F',file_NO);
fid_F = fopen(fname_F,'w');

fname_P = sprintf('%s_%d','./data/P',file_NO);
fid_P = fopen(fname_P,'w');

fname_I = sprintf('%s_%d','./data/I',file_NO);
fid_I = fopen(fname_I,'w');

fname_O = sprintf('%s_%d','./data/SO',file_NO);
fid_O = fopen(fname_O,'w');

fname_N = sprintf('%s_%d','./data/SN',file_NO);
fid_N = fopen(fname_N,'w');

fname_A = sprintf('%s_%d','./data/SA',file_NO);
fid_A = fopen(fname_A,'w');

fname_G = sprintf('%s_%d','./data/SG',file_NO);
fid_G = fopen(fname_G,'w');

if Ind_B == 1
    fname_B = sprintf('%s_%d','./data/TB',file_NO);
    fid_B = fopen(fname_B,'w');
end

if Ind_T == 1
    fname_T = sprintf('%s_%d','./data/TT',file_NO);
    fid_T = fopen(fname_T,'w');
end

fname_tf = sprintf('%s_%d','./data/theta_f',file_NO);
fid_tf = fopen(fname_tf,'w');

fname_tp = sprintf('%s_%d','./data/theta_p',file_NO);
fid_tp = fopen(fname_tp,'w');

for iw = 1:imax
    fprintf(fid_F, '%12.10e ', F2D(iw,:));
    fprintf(fid_F, '\n');
    
    fprintf(fid_P, '%12.10e ', P2D(iw,:));
    fprintf(fid_P, '\n');
    
    fprintf(fid_I, '%12.10e ', I2D(iw,:));
    fprintf(fid_I, '\n');
    
    fprintf(fid_O, '%12.10e ', SO2D(iw,:));
    fprintf(fid_O, '\n');
    
    fprintf(fid_N, '%12.10e ', SN2D(iw,:));
    fprintf(fid_N, '\n');
    
    fprintf(fid_A, '%12.10e ', SA2D(iw,:));
    fprintf(fid_A, '\n');
    
    fprintf(fid_G, '%12.10e ', SG2D(iw,:));
    fprintf(fid_G, '\n');
    
    if Ind_B == 1
        fprintf(fid_B, '%12.10e ', TB2D(iw,:));
        fprintf(fid_B, '\n');
    end
    
    if Ind_T == 1
        fprintf(fid_T, '%12.10e ', TT2D(iw,:));
        fprintf(fid_T, '\n');
    end
        
    fprintf(fid_tf, '%12.10e ', theta_f_2D(iw,:));
    fprintf(fid_tf, '\n');
    
    fprintf(fid_tp, '%12.10e ', theta_p_2D(iw,:));
    fprintf(fid_tp, '\n');
end

fclose(fid_F);
fclose(fid_P);
fclose(fid_I);
fclose(fid_O);
fclose(fid_N);
fclose(fid_A);
fclose(fid_G);

if Ind_B == 1
    fclose(fid_B);
end

if Ind_T == 1
    fclose(fid_T);
end

fclose(fid_tf);
fclose(fid_tp);

tmp = 0;

end
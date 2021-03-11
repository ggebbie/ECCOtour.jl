addpath /home/swijffels/toolbox/susan
addpath /home/swijffels/toolbox/seawater
addpath ~/toolbox/csirolib
addpath ~/work/argo/matlab

cd /home/swijffels/work/argo/gridNSF/ecco4/

%% interpolates S, T and PT, and velocity onto or sigma1 analysis grid - and gets transports between layers
% set the grid here - same as make_sigsurf_dedrift_2019.m
 
sig1grid =[24:0.05:31, 31.02:0.02:33];
sig1grid = sig1grid(1:3:length(sig1grid));
sig1grid = sig1grid(:);

% tname = ['/data/ecco-v4r3-interp/monthly/THETA.2011.nc'];
% use WHOI generated versions instead
tname = ['/data/ecco-v4r4/interp_monthly/THETA/2000/THETA_2000_01.nc'];
nyears = 2017-1992 +1;  
 nx = 720;ny=360; nz = 50; nt = nyears*12;  % year range

% 
% yi = -70:1:70;   % Argo grid
% xi = 0:2:360;
% xi=xi(:);yi=yi(:);

 cp = sw_cp(34.6,15,0);
 rho = sw_dens(34.6,15,0);
 ppi = 0:.25:2000; ppi = ppi(:); % transport integration grid
%% keep grid in Ecco files
its = 0;
dnum = [];
for iyear = 1999:2017 %1993:2017;
    cyear = int2str(iyear)
    [si,pri,pti,vi,ui]=deal(NaN(nx,ny,length(sig1grid),12));   % cannot do all years  in one go! Not enough memory on atalanta: do in 12 months blocks
    [trv,trs,trh]=deal((1+1i)*NaN(nx,ny,length(sig1grid),12));   % cannot do all years  in one go! Not enough memory on atalanta: do in 12 months blocks
    
    if its == 0 
        x = ncread(tname,'longitude');
        x = x(:);
        
        
        xs = rem(x+360,360);
        [xs,isx]=sort(xs);
        
%         ibad = abs(xs - 180) < .5;
%         isx(ibad) = [];  % bad data interpolation near dateline
%         xs(ibad) = [];
%          xs = [-.25;xs;360.25];
%          isx = [find(x==-.25);isx;find(x==.25)];
        
        y = ncread(tname,'latitude');
        y  = y(:);
%        z = ncread(tname,'dep');
        z = -1*ncread(tname,'Z');  % in meters
        xi=xs;yi = y;
    end
 
%    
%     pt1 = ncread(tname,'THETA');
%     s1 = ncread(sname,'SALT');
%     
 
    
    for it = 1:12
        cmon = int2str(it+100);
        cmon(1) = [];
     tname = ['/data/ecco-v4r4/interp_monthly/THETA/',cyear,'/THETA_',cyear,'_',cmon,'.nc']
     sname = ['/data/ecco-v4r4/interp_monthly/SALT/',cyear,'/SALT_',cyear,'_',cmon,'.nc']; 
     uname = ['/data/ecco-v4r4/interp_monthly/EVEL/',cyear,'/EVEL_',cyear,'_',cmon,'.nc']; 
     vname = ['/data/ecco-v4r4/interp_monthly/NVEL/',cyear,'/NVEL_',cyear,'_',cmon,'.nc']; 
     mday = datenum(iyear,it,15,12,0,0); % middle of the month
        dnum = [dnum;mday];  
     
     % read in data to temp fields
      [s2,pt2,u2,v2]=deal(NaN(length(isx),length(yi),length(z)));
         pt1 = ncread(tname,'THETA');
         s1 = ncread(sname,'SALT'); 
         u1 = ncread(uname,'EVEL');  % m/s
         v1 = ncread(vname,'NVEL'); 
         imask = s1==0;
         pt1(imask)=NaN;
         s1(imask)=NaN;
         u1(imask)=NaN;
         v1(imask)=NaN;
          s2  =  s1(isx,:,:);
          pt2  =  pt1(isx,:,:);
          u2  =  u1(isx,:,:);
          v2  =  v1(isx,:,:);
          
 %        for iz = 1:length(z)
%             s2(:,:,iz) = interp2(xs,y,squeeze(s1(isx,:,iz))',xi,yi')';
%             pt2(:,:,iz) = interp2(xs,y,squeeze(pt1(isx,:,iz))',xi,yi')';
%             u2(:,:,iz) = interp2(xs,y,squeeze(u1(isx,:,iz))',xi,yi')';
%             v2(:,:,iz) = interp2(xs,y,squeeze(v1(isx,:,iz))',xi,yi')';
%           s2(:,:,iz) =  s1(isx,:,iz);
%           pt2(:,:,iz) =  pt1(isx,:,iz);
%           u2(:,:,iz) =  u1(isx,:,iz);
%           v2(:,:,iz) =  v1(isx,:,iz);
%          
 %       end
        % calculate sigma and interpolate pt and S onto the levels we want
        for iy = 1:length(yi)
                pp=sw_pres(z,yi(iy));
                pp(1) = -1; 
                for ix = 1:length(xi)
                
                sig1 = sw_pden(squeeze(s2(ix,iy,:)),squeeze(pt2(ix,iy,:)),0*ones(size(z)),1000)-1000;
                ss = squeeze(s2(ix,iy,:));
                tt = squeeze(pt2(ix,iy,:));
                uu = squeeze(u2(ix,iy,:));
                vv = squeeze(v2(ix,iy,:));

               ig = ~isnan(ss);
                if sum(ig) > 3
                     
                    sig11 = sig1(ig) + [1:sum(ig)]'*2e-6;
%                     si(ix,iy,:,it+its) = interp1q(sig11,ss(ig),sig1grid);
%                     pti(ix,iy,:,it+its) = interp1q(sig11,tt(ig),sig1grid);
%                     pp = sw_pres(z,yi(iy));  % convert to pressure
%                     pri(ix,iy,:,it+its) = interp1q(sig11,pp(ig),sig1grid);
%                     ui(ix,iy,:,it+its) = interp1q(sig11,uu(ig),sig1grid);
%                     vi(ix,iy,:,it+its) = interp1q(sig11,vv(ig),sig1grid);
                    
                    si(ix,iy,:,it) = interp1q(sig11,ss(ig),sig1grid);
                    pti(ix,iy,:,it) = interp1q(sig11,tt(ig),sig1grid);
                    pp = sw_pres(z,yi(iy));  % convert to pressure
                    pri(ix,iy,:,it) = interp1q(sig11,pp(ig),sig1grid);
                    ui(ix,iy,:,it) = interp1q(sig11,uu(ig),sig1grid);
                    vi(ix,iy,:,it) = interp1q(sig11,vv(ig),sig1grid);   
                    
                    % do tranports:
                    ssi = interp1q(pp,ss,ppi);
                    ptti = interp1q(pp,tt,ppi);
                    vvi = interp1q(pp,vv,ppi);
                    uui = interp1q(pp,uu,ppi);
                     vtri = uui + 1i*vvi;
                     stri = ssi.*vtri*rho;
                     htri = ptti.*vtri*rho*cp;
                     
                    icount = 0;
                    prsig = squeeze(pri(ix,iy,:,it));
                    igsig = find(~isnan(prsig))';
                    if ~isempty(igsig) 
                    for is = igsig
                        
                        if icount == 0
                            pupp = 0;  % surface for top layer
                        else
                            pupp = prsig(is-1);
                        end
                        icount = icount+1;
                        
                        ii = ppi <  prsig(is) & ppi >= pupp;
                        if sum(ii)
                            %   plot(vtr2(ii),p2(ii),'g.',squeeze(climv.ug(ixx,iyy,:)),pgrid),axis ij,pause
                            trv(ix,iy,is,it) = nansum(vtri(ii))*.25;
                            trs(ix,iy,is,it) = nansum(stri(ii))*.25;
                            trh(ix,iy,is,it) = nansum(htri(ii))*.25;
                        end
                    end
                    end
%                   plot(cumsum(repnan(vtri,0))*.25,ppi,cumsum(squeeze(repnan(trv(ix,iy,:,it),0))),prsig,'o')
%                   axis ij
%                   pause
                end
                
                end
             [num2str(yi(iy)),' ',datestr(mday)]
        end
        
    end
 %   its = its+12
  %datestr(dnum(end))
  yrgrid = iyear + [0.5:1:(12)]/12;

readme = 'made from JPL generated eccoV4r4  by SandVELonSigma1_eccov4r4_fine_Sep2020.m  ';

sname = ['/home/swijffels/work/argo/gridNSF/ecco4/gridonSig1_ecco4r4_trans_',cyear,'.mat']
save -v7.3 junk.mat xi yi yrgrid si pri pti ui vi sig1grid dnum yrgrid readme trv trs trh
eval(['!mv junk.mat ',sname]);

end
  

exit
 
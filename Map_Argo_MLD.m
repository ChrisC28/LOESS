list_var={'depth','pdens','temp','psal','drhodz_base','NT_base','NS_base'}

for ivar =1:length(list_var)
  var=list_var{ivar};
  disp(['work on: ' var])
  x=0:.5:360;
  y=-90:.5:90;
  % nb: 12 car chaque mois
  eval(['MLfitted.' var '=NaN*ones(12,length(x),length(y));'])
  eval(['MLfitted.' var '_Nbpts=NaN*ones(12,length(x),length(y));'])
  eval(['MLfitted.' var '_RadiusFit=NaN*ones(12,length(x),length(y));'])
  for i=1:12
    disp(['Mois:' int2str(i)])

    % definition de la grille
    [xx , yy]=meshgrid(x,y);
    xg=reshape(xx,length(x)*length(y),1);
    yg=reshape(yy,length(x)*length(y),1);

    % Selection des donnÃ©es:
    iok=find(Month==i);  
    xd=ML.lon(iok); 
    yd=ML.lat(iok);
    tims=ML.date(iok)-datenum(1900,1,1);
    eval(['dd=ML.' var  '(iok);']);

    % on enleve les mauvaises valeurs
    inotok=find(isnan(xd+yd+dd+tims) | xd<-180 | yd>90 | yd<-90 | xd>180);
    dd(inotok)=[];
    xd(inotok)=[];
    yd(inotok)=[];
    tims(inotok)=[];
    xd(find(xd<0))=xd(find(xd<0))+360;

    % Fitting: 
    dw=[];
    fns=[0 0 0 0];
    bar=0;
    scrn=0;
    fopt = [1 2 5 11 41];
    fval = {300,2500,1,1,0};
    [mn,nq,rq] = loess_map(xg,yg,xd,yd,dd,tims,dw,fns,bar,scrn,fopt,fval);
    Result=reshape(mn,length(y),length(x));
    MaxR=reshape(rq,length(y),length(x));
    Nbpt=reshape(nq,length(y),length(x));
    eval(['MLfitted.' var '(i,:,:)=Result'';'])
    eval(['MLfitted.' var '_Nbpts(i,:,:)=Nbpt'';'])
    eval(['MLfitted.' var '_RadiusFit(i,:,:)=MaxR'';'])
    save ../Matrix/MLfitted.mat MLfitted
  end
  save ../Matrix/MLfitted.mat MLfitted
end
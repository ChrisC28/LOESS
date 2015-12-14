clear all;
rand('seed',100);

nPoints = 2000;
%0:.5:360
x_west  =  0.0;
x_east  = +360.0;
y_south =  00.0;
y_north = +90.0;

x_points = x_west  + (x_east-x_west)  .* rand([nPoints,1]);
y_points = y_south + (y_north-y_south).* rand([nPoints,1]);

sigma_x = (x_east  - x_west)  * 0.25
sigma_y = (y_north - y_south) * 0.25
z_points = (x_points - 0.5*(x_west+x_east)).*exp(-( ( x_points - 0.5*(x_west+x_east)   ).^2)/(2.0*sigma_x*sigma_x)  ...
		                                 -( ( y_points - 0.5*(y_south+y_north) ).^2)/(2.0*sigma_y*sigma_y) );
size(z_points)

x = x_west:0.5:x_east;
y = y_south:0.5:y_north;

[x_grid,y_grid] = meshgrid(x,y);


z_grid = (x_grid - 0.5*(x_west+x_east)).*exp(-( ( x_grid - 0.5*(x_west+x_east)   ).^2)/(2.0*sigma_x*sigma_x)  ...
		                             -( ( y_grid - 0.5*(y_south+y_north) ).^2)/(2.0*sigma_y*sigma_y) );

x_grid=reshape(x_grid,length(x)*length(y),1);
y_grid=reshape(y_grid,length(x)*length(y),1);


dw=[];
fns=[0 0 0 0];
bar=0;
scrn=0;
fopt = [1 2,5,11,13];
fval = {100,500,0,2,0};


times = 1.0 * ones(size(z_points));

[mn,nq,rq] = loess_map(x_grid,y_grid,x_points,y_points,z_points,times,dw,fns,bar,scrn,fopt,fval);    

Result=reshape(mn,length(y),length(x));
MaxR=reshape(rq,length(y),length(x));
Nbpt=reshape(nq,length(y),length(x));
    
x_grid=reshape(x_grid,length(y),length(x));
y_grid=reshape(y_grid,length(y),length(x));

size(x_grid)
size(y_grid)
size(z_grid)

figure(1)
contourf(x_grid,y_grid,z_grid)
hold on
plot(x_points,y_points,'r*')
figure(2)
contourf(x_grid,y_grid,Result)





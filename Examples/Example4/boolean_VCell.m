period = 1;
stim_duration = 0.1;
duration = 100;
stim_s = 0:period:duration;
stim_e = (0:period:duration) + stim_duration;

origin = [1.82895,-3.20716,-12.4074];
size_xyz = [30.9448,40.9563,28.1864];
mesh_xyz = [60,79,55];
%mesh_xyz = [10,10,10];
desired_sampling = mesh_xyz * 1;
activation = zeros(desired_sampling);
%activation(4:6,4:6,4:6) = 1;
activation(25:35,35:45,22:31) = rand(11,11,10);

activation = activation(:);

x = linspace(origin(1),size_xyz(1),desired_sampling(1));
y = linspace(origin(2),size_xyz(2),desired_sampling(2));
z = linspace(origin(3),size_xyz(3),desired_sampling(3));
[X,Y,Z] = meshgrid(x,y,z);
X = X(:);
Y = Y(:);
Z = Z(:);
step_x = (x(2)-x(1))/2;
step_y = (y(2)-y(1))/2;
step_z = (y(2)-y(1))/2;
all = [];
for i = 1:size(X)
    if activation(i) < 1e-3
        continue
    end
    xc = ['(x>',num2str(X(i)-step_x),')&&(x<=',num2str(X(i)+step_x),')'];
    yc = ['(y>',num2str(Y(i)-step_y),')&&(y<=',num2str(Y(i)+step_y),')'];
    zc = ['(z>',num2str(Z(i)-step_z),')&&(z<=',num2str(Z(i)+step_z),')'];
    cc = ['(',num2str(activation(i)),'*(',xc,'&&',yc,'&&',zc,'))'];
    if isempty(all)
        all = [cc];
    else
        all = [all,'+',cc];
    end
end

time_logic = [];
for i = 1:size(stim_s,2)
    tc = ['(t>=',num2str(stim_s(i)),')&&(t<',num2str(stim_e(i)),')'];
    if isempty(time_logic)
        time_logic = [tc];
    else
        time_logic = [time_logic,'+',tc];
    end
end

fileID = fopen('test.txt','w');
fprintf(fileID,'%s\n', all)
fprintf('%s\n', time_logic)

%% TIRF

on_raster = zeros(21,21);
on_raster(6:10:end,6:10:end) = 1;
d =  @(lambda,n1,n2,theta) lambda/4/pi * (n2^2 * sind(theta)^2 - n1^2)^(-1/2);
wave_z = @(lambda,n1,n2,theta,z) exp(-z/d(lambda,n1,n2,theta));
lambda = .445;
n1 = 1.33;
n2 = 1.52;
n = 1.33;
theta = 79;
w0 = sqrt(2) * 0.325 * lambda / (sqrt(2) * 1.4^0.91);
zr = pi*w0^2*n/lambda;
w = @(z) w0*sqrt(1 + (z/zr).^2);
profile = @(r,z) exp(-2*r.^2./(w(z).^2));

d = lambda/4/pi * (n2^2 * sind(theta)^2 - n1^2)^(-1/2);
w0 = sqrt(2) * 0.325 * lambda / (sqrt(2) * 1.4^0.91);
zr = pi*w0^2*n/lambda;
x0b = 30;
y0b = 30;
z0 = 0.27972204;


rx = -0.95:0.1:1;
ry = -0.95:0.1:1;
%rx = -0.95:0.5:1;
%ry = -0.95:0.5:1;
[RX,RY] = meshgrid(rx,ry);
RX = RX(:);
RY = RY(:);

eq_vect = [];
for i = 1:size(RX,1)
    x0 = x0b + RX(i);
    y0 = y0b + RY(i);
    s_o = ['exp(-(z-(',num2str(z0),'))/',num2str(d),')*',...
        'exp(-2*((x-(',num2str(x0),'))*(x-(',num2str(x0),'))+(y-(',...
        num2str(y0),'))*(y-(',num2str(y0),')))/',...
        num2str(w0^2),'/(1+(z-(',num2str(z0),'))*(z-(',...
        num2str(z0),'))/',num2str(zr^2),'))'];
    if isempty(eq_vect)
        eq_vect = s_o
    else
        eq_vect = [eq_vect,'+',s_o];
    end
end
fileID = fopen('TIRF.txt','w');
%fprintf(fileID,'%s\n', all)
fprintf(fileID,'%s\n', ['(t<=0.1)*(',eq_vect],')')

%% 1P

on_raster = zeros(21,21);
on_raster(6:10:end,6:10:end) = 1;
%1P specific stuff%
lambda = .450;
n=1.33;
w0 = sqrt(2) * 0.325 * lambda / (sqrt(2) * 1.4^0.91);
zr = pi*w0^2*n/lambda;
w = @(z) w0*sqrt(1 + (z/zr).^2);
I = @(r,z) (w0./w(z)).^2 .* exp(-2*r.^2./(w(z).^2));

x0b = 30;
y0b = 30;
z0 = 0.27972204;



rx = -0.95:.1:1;
ry = -0.95:.1:1;
%rx = -0.95:0.5:1;
%ry = -0.95:0.5:1;
[RX,RY] = meshgrid(rx,ry);
RX = RX(:);
RY = RY(:);

eq_vect = [];
for i = 1:size(RX,1)
    x0 = x0b + RX(i);
    y0 = y0b + RY(i);
    s_o = ['(1/(1+(z-(',num2str(z0),'))*(z-(',...
        num2str(z0),'))/',num2str(zr^2),'))*(exp(-2*((x-(',num2str(x0),'))*(x-(',num2str(x0),'))+(y-(',...
        num2str(y0),'))*(y-(',num2str(y0),')))/',...
        num2str(w0^2),'/(1+(z-(',num2str(z0),'))*(z-(',...
        num2str(z0),'))/',num2str(zr^2),')))'];
    if isempty(eq_vect)
        eq_vect = s_o;
    else
        eq_vect = [eq_vect,'+',s_o];
    end
end
fileID = fopen('1P_offset.txt','w');
fprintf(fileID,'%s\n', ['(t<=0.1)*(',eq_vect,')'])

%% 2P

on_raster = zeros(21,21);
on_raster(6:10:end,6:10:end) = 1;
%1P specific stuff%
lambda = .900;
n=1.33;
w0 = sqrt(2) * 0.325 * lambda / (sqrt(2) * 1.4^0.91);
zr = pi*w0^2*n/lambda;
w = @(z) w0*sqrt(1 + (z/zr).^2);
I = @(r,z) (w0./w(z)).^2 .* exp(-2*r.^2./(w(z).^2));

x0b = 30;
y0b = 30;
z0 = 0.27972204;



rx = -0.95:.1:1;
ry = -0.95:.1:1;
%rx = -0.95:0.5:1;
%ry = -0.95:0.5:1;
[RX,RY] = meshgrid(rx,ry);
RX = RX(:);
RY = RY(:);

eq_vect = [];
for i = 1:size(RX,1)
    x0 = x0b + RX(i);
    y0 = y0b + RY(i);
    s_o = ['(1/(1+(z-(',num2str(z0),'))*(z-(',...
        num2str(z0),'))/',num2str(zr^2),'))^2*(exp(-2*((x-(',num2str(x0),'))*(x-(',num2str(x0),'))+(y-(',...
        num2str(y0),'))*(y-(',num2str(y0),')))/',...
        num2str(w0^2),'/(1+(z-(',num2str(z0),'))*(z-(',...
        num2str(z0),'))/',num2str(zr^2),'))^2)'];
    if isempty(eq_vect)
        eq_vect = s_o;
    else
        eq_vect = [eq_vect,'+',s_o];
    end
end
fileID = fopen('2P_offset.txt','w');
fprintf(fileID,'%s\n', ['(t<=0.1)*(',eq_vect,')'])
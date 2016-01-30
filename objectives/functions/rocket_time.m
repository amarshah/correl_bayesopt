function[ output ] = rocket_time( x )

ff = 0.95+0.03*x(1); % fuel fraction
an = x(2)*90; % angle of burn
y0 = 1+0.02*x(3);

%rocket.m
%program to solve the variable mass rocket equation in the
%presence of a gravitational force from a massive body
%clear;
G=6.67e-11;                    %universal gravitational constant (Nm^2/kg^2)
R=6.37e6;                      %unit of distance (m) - massive body radius
m0=5.98e24;                    %unit of mass (kg) - could be the massive body
M=5.98e24/m0;                  %massive body mass in units of m0
tau=sqrt(R^3/(G*m0));          %unit of time in seconds
v0=R/tau; F0=m0*R/tau^2;       %unit of speed and unit of force
mi=2.8e6/m0;                   %initial payload+fuel mass in units of m0
mp=mi*ff;             %fuel fraction, and fuel mass
Thrust=1.5*(G*mi*M/R^2)*m0^2;  %Let the Thrust be # times initial mass weight
Thrust=Thrust/F0;              %Thrust in units of force
u=0.35;                        %gas exhaust velocity in units of v0
th=an*2*pi/360;         %angle of burn determines launch angle
ux=u*cos(th); uy=u*sin(th);    %exhaust velocity components in units of v0
alpha=(Thrust/u);              %alpha in units of m0/tau
mf=mi-mp;                      %final mass is payload mass (after fuel burnout)
tfmax=(mi-mf)/alpha;           %fuel burnout time in units of tau
tmax=1.5e3;                 %simulation run time
x0=0;vx0=0;vy0=0;         %initial positions, speeds
ic1=[x0;vx0;y0;vy0;mi];        %initial conditions: position, velocity, rocket mass
%Use MATLAB's Runge-Kutta (4,5) formula (uncomment/comment as needed)
%opt=odeset('AbsTol',1.e-8,'RelTol',1.e-5);%user set Tolerances
%[t,w]=ode45('rocket_der',[0.0,tmax],ic1,opt,alpha,ux,uy,M,tfmax);%with set tolerance
[t,w]=ode45('rocket_der',[0.0,tmax],ic1,[],alpha,ux,uy,M,tfmax);%default tolerance
L=2.5*sqrt(x0^2+y0^2);         %window size
h=[0:0.025:2*pi];
x=cos(h);y=sin(h);             %massive body perimeter
%plot(x,y,'g'); hold on        %massive body plot if needed here
%plot(w(:,1),w(:,3))           %use this to plot all the x,y points calculated
n=length(t);                   %size of the time array
angle = 0;
for i=1:n                      %Loop to pick the points that lie above ground
  if sqrt(w(i,1)^2+w(i,3)^2) >= 0.99 %ground is 1.0 R, include points slightly below
     nn=i;
     t1(i)=t(i);
     x1(i)=w(i,1); 
     y1(i)=w(i,3);
     vx1(i)=w(i,2); 
     vy1(i)=w(i,4);
     
     if i>1
         angle_iminus1 = atan(y1(i-1)/x1(i-1));
         angle_i       = atan(y1(i)/x1(i));
         angle = angle + abs(angle_i-angle_iminus1);
     end
  else
      break;
  end
end

time = 2*t1(nn)/tmax;
angle = min(2,angle /2/pi);
if time>=3
    time=0;
    angle=0;
end

time = time + 3;

launch_angle_penalty = (1-an/90);

output = [time, angle, launch_angle_penalty]';

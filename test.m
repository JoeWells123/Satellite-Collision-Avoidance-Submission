initialTime = datetime(2022,12,14,1,4,0);
initialTimeJD = juliandate(initialTime);

a1 = 7000e3;
e1 = 0;
i1 = 45;
raan1 = 90;
w1 = 30;
theta1 = 0;
mu = 398600.4418e9;

t1 = 2*pi*sqrt((a1^3)/mu);
t1 = 0;
a2 = 10000e3;
aTransfer = (a1 + a2)/2;
a2 = 10000e3;

aTransfer = (a1 + a2)/2;
tTransfer = pi*sqrt((aTransfer^3)/mu);
t2 = t1 + tTransfer;
tf = t2;
disp(tf)
v1 = sqrt(mu/a1);
vTransfer1 = sqrt(2*mu*((1/a1) - (1/(2*aTransfer))));
deltav1 = vTransfer1 - v1;
v2 = sqrt(mu/a2);
vTransfer2 = sqrt(2*mu*((1/a2) - (1/(2*aTransfer))));
deltav2 = v2 - vTransfer2;

thetaStar = 0;

deltav1Direction = [- cosd(raan1)*sind(w1 + thetaStar) - sind(raan1)*cosd(i1)*cosd(w1 + thetaStar); ...
    - sind(raan1)*sind(w1 + thetaStar) + cosd(raan1)*cosd(i1)*cosd(w1 + thetaStar); ...
    sind(i1)*cosd(w1 + thetaStar)];
thetaStar = 180;
deltav2Direction = [- cosd(raan1)*sind(w1 + thetaStar) - sind(raan1)*cosd(i1)*cosd(w1 + thetaStar); ...
    - sind(raan1)*sind(w1 + thetaStar) + cosd(raan1)*cosd(i1)*cosd(w1 + thetaStar); ...
    sind(i1)*cosd(w1 + thetaStar)];


g0 = 9.81;
Isp = 400;
ve = g0*Isp;
mDot = 500;
m0 = 1000;

m1 = m0/exp(deltav1/ve);
burnDuration1 = (m0-m1)/mDot;
m2 = m1/exp(deltav2/ve);
burnDuration2 = (m1-m2)/mDot;


model = 'hohmannTransfer';
open_system(model);

simOut = sim(model);

positionTT = timeseries2timetable(simOut.yout{1}.Values);
attitudeTT = timeseries2timetable(simOut.yout{2}.Values);

startTime = initialTime;
stopTime = startTime + seconds(tf);
sampleTime = 60;
sc = satelliteScenario(startTime,stopTime,sampleTime);


sat = satellite(sc,positionTT,Name = "Spacecraft");

pointAt(sat,attitudeTT);

% Launch a satellite scenario viewer
v = satelliteScenarioViewer(sc,CameraReferenceFrame = "Inertial");
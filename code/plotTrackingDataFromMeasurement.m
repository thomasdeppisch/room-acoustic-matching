% This script reads and displays tracking information obtained from
% OptiTrack.
% td, 2025

clear all
close all

addpath("lib/")

dataStruct = load('../data/rooms/hallway_2p5m_femaleSpeech_combinedMove.mat');

%% plot tracking data obtained from a measurement with optitrack
printPos = true;

% prepare plot
plotLimM = 3;

% plot 3 different views of the same 3d plot
hPosMarker = {};
hRotArrow = {};

trackingFig = figure;
subplot(131)
hPosMarker{1} = plot3(0,0,0,'o');
hold on
hRotArrow{1} = quiver3(0,0,0,1,0,0);
grid on
axis([-plotLimM plotLimM -plotLimM plotLimM -plotLimM plotLimM])
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
view([-90 90])
title('top view')

subplot(132)
hPosMarker{2} = plot3(0,0,0,'o');
hold on
hRotArrow{2} = quiver3(0,0,0,1,0,0);
grid on
axis([-plotLimM plotLimM -plotLimM plotLimM -plotLimM plotLimM])
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
view([-90 0])
title('behind view')

subplot(133)
hPosMarker{3} = plot3(0,0,0,'o');
hold on
hRotArrow{3} = quiver3(0,0,0,1,0,0);
grid on
axis([-plotLimM plotLimM -plotLimM plotLimM -plotLimM plotLimM])
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
view([-37.5-90 30])
title('side view')

set(trackingFig,'units','centimeters','position',[5 5 40 10],'paperunits','centimeters','papersize',[40 10])

%% plot tracking
numSamples = size(dataStruct.translationXYZCart,1);
rb = {};
timeStepMs = 100;
sampleStep = round(timeStepMs/1000*dataStruct.fs);

for ii = 1:sampleStep:numSamples
    rb.x = dataStruct.translationXYZCart(ii,1);
    rb.y = dataStruct.translationXYZCart(ii,2);
    rb.z = dataStruct.translationXYZCart(ii,3);
    rb.yawRad = dataStruct.rotationYawPitchRollRad(ii,1);
    rb.pitchRad = dataStruct.rotationYawPitchRollRad(ii,2);
    rb.rollRad = dataStruct.rotationYawPitchRollRad(ii,3);
    plotRigidBody(hPosMarker,hRotArrow,rb,printPos)
    pause(0.01);
end


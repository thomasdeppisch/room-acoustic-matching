function plotRigidBody(hPosMarker,hRotArrow,rigidBody,printPos)

x = rigidBody.x;
y = rigidBody.y;
z = rigidBody.z;

yawRad = rigidBody.yawRad;
pitchRad = rigidBody.pitchRad;
rollRad = rigidBody.rollRad;

if printPos
    fprintf( 'X:%0.1f mm  ', x * 1000 )
    fprintf( 'Y:%0.1f mm  ', y * 1000 )
    fprintf( 'Z:%0.1f mm\n', z * 1000 )
    fprintf( 'eulerX:%0.1f deg  ', rollRad*180/pi )
    fprintf( 'eulerY:%0.1f deg  ', pitchRad*180/pi )
    fprintf( 'yaw:%0.1f deg\n', yawRad*180/pi )
end

dir = [1 0 0 1];
rotMtx = makehgtform('xrotate',rollRad,'yrotate',pitchRad,'zrotate',yawRad);
rotDirCart = rotMtx * dir';

numSubplots = length(hPosMarker);

% update all subplots
for ii = 1:numSubplots
    hPosMarker{ii}.XData = x;
    hPosMarker{ii}.YData = y;
    hPosMarker{ii}.ZData = z;
    
    hRotArrow{ii}.XData = x;
    hRotArrow{ii}.YData = y;
    hRotArrow{ii}.ZData = z;
    
    hRotArrow{ii}.UData = rotDirCart(1);
    hRotArrow{ii}.VData = rotDirCart(2);
    hRotArrow{ii}.WData = rotDirCart(3);
end

drawnow
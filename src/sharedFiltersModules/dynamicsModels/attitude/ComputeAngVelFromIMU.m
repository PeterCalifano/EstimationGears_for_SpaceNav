function [dEstAngVel] = ComputeAngVelFromIMU(dMeasAngVel, dGyrosBias, dDCM_IMUfromNavFrame, dNavFrameAngVel)
arguments
    dMeasAngVel          (3,1) double {isvector, isnumeric}
    dGyrosBias           (3,1) double {isvector, isnumeric} 
    dDCM_IMUfromNavFrame (3,3) double {ismatrix, isnumeric} = zeros(3,3)
    dNavFrameAngVel      (3,1) double {isvector, isnumeric} = zeros(3,1)
end

dNonInertialAngVel = coder.nullcopy(zeros(3,1));

if any(dNavFrameAngVel > 0, "all")
    assert(all(dNavFrameAngVel ~= 0, 'all'));
end

% Compute angular velocity removing bias
dEstAngVel = coder.nullcopy(zeros(3,1));
dEstAngVel(1:3) = dMeasAngVel(1:3) - dGyrosBias(1:3);

% Remove non inertial term (if not zero)
if any(dDCM_IMUfromNavFrame > 0, "all") && any(dNavFrameAngVel, 'all')
    dNonInertialAngVel(1:3) = dDCM_IMUfromNavFrame * dNavFrameAngVel;
    dEstAngVel = dEstAngVel - dNonInertialAngVel;
end

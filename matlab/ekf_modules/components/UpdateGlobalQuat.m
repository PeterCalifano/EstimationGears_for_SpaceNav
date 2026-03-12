function [dUpdatedQuat] = UpdateGlobalQuat(dQuat, dDeltaErrQuat) %#codegen
arguments
    dQuat       (3,1) double 
    dDeltaErrQuat (4,1) double
end

dUpdatedQuat = coder.nullcopy(zeros(4,1));

% Extract components
dqw1 = dQuat(1); 
dqx1 = dQuat(2); 
dqy1 = dQuat(3); 
dqz1 = dQuat(4);

dqw2 = dDeltaErrQuat(1); 
dqx2 = dDeltaErrQuat(2); 
dqy2 = dDeltaErrQuat(3); 
dqz2 = dDeltaErrQuat(4);

% Quaternion multiplication (dQuat * dDeltaErrQuat)
dqw_new = dqw1*dqw2 - dqx1*dqx2 - dqy1*dqy2 - dqz1*dqz2;
dqx_new = dqw1*dqx2 + dqx1*dqw2 + dqy1*dqz2 - dqz1*dqy2;
dqy_new = dqw1*dqy2 - dqx1*dqz2 + dqy1*dqw2 + dqz1*dqx2;
dqz_new = dqw1*dqz2 + dqx1*dqy2 - dqy1*dqx2 + dqz1*dqw2;

% Concat updated quaternion
dUpdatedQuat(:) = [dqw_new; dqx_new; dqy_new; dqz_new];
% dUpdatedQuat = dUpdatedQuat / norm(dUpdatedQuat);

end

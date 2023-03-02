function K = calibrationParameterstoMatrix(calibration)
% convert calibration parameters to something useful and consistent
if isa(calibration,'cameraParameters') % if we have cameraparameters from calibration toolbox
    K.inv = calibration.IntrinsicMatrix'; % the intrinsic paramters are given as the inverse of the transpose
    K.fwd = inv(K.inv); % invert them for the forward
elseif isstruct(calibration); % otherwise throw an error
    K = calibration;
elseif ismatrix(calibration) % if we are given just a matrix
    K.fwd = calibration; % assume it's the forward transform
    K.inv = inv(calibration); % and invert it for the reverse transform
else
    error('unknown calibration input for conversion')
end
    
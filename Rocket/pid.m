Thrust = 20;
moi = 0.09;
TVCTheta = [0,0];

horizontalForce = Thrust * sin(TVCTheta);
angularAcceleration = horizontalForce / moi;
theta = integrate2(angularAcceleration);

blk = tunablePID('pdblock','PID');
blk.Kp.Value = 4;        % initialize Kp to 4
blk.Ki.Value = 0.1;
blk.Kd.Value = 0.7;      % initialize Kd to 0.7
blk.Tf.Value = 0.01;     % set parameter Tf to 0.01
blk.Tf.Free = false;     % fix parameter Tf to this value
blk;

TVCTheta = blk;
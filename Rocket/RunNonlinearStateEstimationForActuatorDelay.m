function final_state = RunNonlinearStateEstimationForActuatorDelay(time, state, Actuator)
  model = 'StateEstimatorForActuatorDelay';
  %cd 'C:\Users\rag\Documents\TVC\Rocket'
  %load_system(model)

  in = Simulink.SimulationInput(model);
  in = in.setVariable('initial_v', state(11:13));
  in = in.setVariable('initial_x', state(14:16));
  in = in.setVariable('initial_w', state(1:3));
  in = in.setVariable('initial_q', state(4:7));
  in = in.setVariable('initial_qi', state(8:10));
  in = in.setVariable('t0', time);
  
  in = in.setVariable('Moments', Actuator);

  simOut = sim(in);
  
  Quaternion_integral = simOut.xFinal{1}.Values;
  Quaternion = simOut.xFinal{3}.Values;
  AngularVelocity = simOut.xFinal{4}.Values;
  temp = simOut.xFinal{2}.Values;
  Xe = temp.Data(end, 1:3);
  Ve = temp.Data(end, 4:6);
  
  final_state = [AngularVelocity.Data(end, :), Quaternion.Data(end, :), Quaternion_integral.Data(end, :), Ve, Xe];
end
# Model of the Mamalian ET cell

Author: Ryan Viertel

## Usage:

data = ET(input,sampling_rate);

input: input vector, if no input then just a vector of zeros

sampling_rate: rate at which the input vector should be sampled. 1000 for milisecond (default)

returns the following struct:

data.T - time vector
data.X - ODE variables at each time step
* data.X(:,1) - Membrane Potential
* data.X(:,2) - nK
* data.X(:,3) - hNaP
* data.X(:,4) - hH
* data.X(:,5) - mLVA
* data.X(:,6) - hLVA
* data.X(:,7) - wBK
* data.X(:,8) - Calcium
* data.X(:,9) - nNew

data.events - list of spike events
data.which - event type
* 1 - spike
* 2 - burst start
 *3 - burst end

data.current - system currents
* data.current(:,1) = transient sodium
* data.current(:,2) = fast potassium
* data.current(:,3) = leak
* data.current(:,4) = persistent sodium
* data.current(:,5) = hyperpolarization activated
* data.current(:,6) = LVA calcium
* data.current(:,7) = HVA calcium
* data.current(:,8) = large conductance potassium

## example

### create the input vector
input = zeros(1,5000);
### run the model
data = ET(input,1000); (or simply "data = ET(input);")
### plot the voltage trace
plot(data.T,data.X(:,1)

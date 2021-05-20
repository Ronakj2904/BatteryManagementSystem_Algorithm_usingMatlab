%%

% function rcValues = tuneModel
%
% This function specifies the resistor and capacitor values that
% you choose to use for either an Rint model or a Thévenin model 
% (see lesson 2.1.3), or an "extended Thévenin model" having 
% additional parallel resistor-capacitor branches.
% 
% If you wish to create an Rint model, simply set rcValues equal
% to the value of R0 (in milliohms).
% If you wish to create a Thévenin model, define rcValues to be 
% a vector having three elements. The first element is R0 (in
% milliohms), the second element is R1 (in milliohms), and the 
% third element is C1 (in kilofarads).
% If you wish to create an extended Thévenin model having 
% additional resistor-capacitor branches, then define rcValues to
% be a vector where the first element is R0 in milliohms, the 
% even-index elements (rcValues(2:2:end)) are resistor values for
% the resistor-capacitor branches, in milliohms, and the remaining
% odd-index elements (rcValues(3:2:end)) are capacitor values for
% the resistor-capacitor branches, in kilofarads.
%
% It is possible to receive full credit for this assignment with
% a well-tuned Thévenin model, but it is also possible to get a 
% much better fit to the data using an extended Thévenin model.
function rcValues = tuneModel

  % BEGIN MODIFYING CODE AFTER THIS
  rcValues = [13.8 ;16; 46];% This is a sample value. You will need to change it.
end  

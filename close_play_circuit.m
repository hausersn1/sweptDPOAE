function close_play_circuit(f1,RP)

%At the end of an experiment, the circuit is removed
invoke(RP,'ClearCOF'); %clear the circuit from the RP2
close(f1); %close the ActX control window
end % of close_play_circuit
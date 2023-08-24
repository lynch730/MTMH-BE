function EN = pulse_lokib(time, charTime, amplitude)
    EN = amplitude.*sqrt(time./charTime).*exp(-time./charTime);
end
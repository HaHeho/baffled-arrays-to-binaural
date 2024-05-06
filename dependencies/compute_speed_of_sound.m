function C = compute_speed_of_sound(temp)

    if isempty(temp)
        C = 343; % in m/s, speed of sound
    else
%         % formlua from Wikipedia
%         C = 331.3 + sqrt(1 + (temp / 273.15));
        
        % formlua with Taylor expansion from Wikipedia
        C = 331.3 + (temp * 0.606);
        
        % I was no able to find a straightforward formula that would also
        % incorporate humidity (it may be not very relevant at all).
    end
end

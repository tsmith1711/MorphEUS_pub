function hextime = get_hex_time()
    % function to get current time and convert to hex for a figure extension
    % use built in now to get current time
    % can't use decimals in hex, so multiply by 10^3 and floor 
    dectime = floor(now*10^3);
    hextime = dec2hex(dectime);
end

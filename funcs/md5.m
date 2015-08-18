function hash = md5(data)
    if ~isa(data,'uint8')
    data1=cast(data,'uint8');
    else
        data1=data;
    end
    x    = java.security.MessageDigest.getInstance('MD5');
    x.update(data1);
    % Hex hash correctly padded. 
    hash = lower(reshape(dec2hex(typecast(x.digest(),'uint8'))',1,[]));
end
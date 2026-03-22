function result = getRMSE( A )
    result = sqrt( sum( sum( A.^2 ) ) / numel( A ) );
end
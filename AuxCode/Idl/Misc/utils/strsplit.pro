function strsplit, stringIn, pattern, _ref_extra=extra

    ON_ERROR, 2  ; return to caller
    RETURN, (n_params() eq 1) ? STRTOK(stringIn, _STRICT_EXTRA=extra) : $
        STRTOK(stringIn, pattern, _STRICT_EXTRA=extra)

end

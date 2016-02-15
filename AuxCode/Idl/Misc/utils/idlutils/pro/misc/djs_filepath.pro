function djs_filepath, filename, root_dir=root_dir, $
 subdirectory=subdirectory, terminal=terminal, tmp=tmp

   if (keyword_set(terminal) OR keyword_set(tmp)) then $
    return, filepath(filename, terminal=terminal, tmp=tmp)

   if (NOT keyword_set(root_dir)) then begin
      if (NOT keyword_set(subdirectory)) then $
       return, filename $
      else $
       return, concat_dir(subdirectory, filename)
   endif else begin
      return, filepath(filename, root_dir=root_dir, $
       subdirectory=subdirectory, terminal=terminal, tmp=tmp)
   endelse

end


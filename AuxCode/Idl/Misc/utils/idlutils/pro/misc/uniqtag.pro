; return elements with unique values of a tag
function uniqtag, str, tag

itag=tag_indx(str[0],tag)
isort=sort(str.(itag))
iuniq=uniq(str[isort].(itag))
return, str[isort[iuniq]]

end

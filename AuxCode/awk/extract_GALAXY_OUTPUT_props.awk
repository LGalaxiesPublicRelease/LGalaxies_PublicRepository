# extract a block starting with a line containing GALAXY_OUTPUT and ending with a }
BEGIN {
flag=0
line=""
d=""
ia=0
b=""
}
/GALAXY_OUTPUT/,/\}/ {   
line=$0
if(flag == 0 && line ~ /\{/) {
	ic=index(line,"{")
	line=substr(line,ic+1)
	flag=1
}
if(line ~ /\}/) {
	ic=index(line,"}")
	if(ic > 0)
	line=substr(line,0,ic-1)
}
if(flag==1 && line != "" && line !~ /^#/ && line !~ /\*\// && line !~ /"GALAXY_OUTPUT"/) 
{
    print line

}
}
# extract a block starting with a line containing GALAXY_OUTPUT and ending with a }
BEGIN {
flag=0
line=""
d=""
}
/GALAXY_OUTPUT/,/\}/ {   
line=$0
if(flag == 0 && line ~ /\{/) {
	ic=index(line,"{")
	line=substr(line,ic+1)
	flag=1
print "GALAXY_OUTPUT {"
}
if(line ~ /\}/) {
	ic=index(line,"}")
	if(ic > 0)
	line=substr(line,0,ic-1)
print "}"
}
if(flag==1 && line != "" && line !~ /^#/) 
{
    isc=index(line,";")
    if(isc == 0)
      d=d" "line
    else {
	    while(isc > 0)
	    {
	      d= d" "substr(line,0,isc-1)
	      # declaration d is now complete, do something with it
#print "<FIELD>"
	      print d
#print "</FIELD>"
	      # reset, see if multiple declarations on one line
	      d=""
	      line=substr(line,isc+1)
	      isc=index(line,";")
	    }
	    d=line
	}
}
}
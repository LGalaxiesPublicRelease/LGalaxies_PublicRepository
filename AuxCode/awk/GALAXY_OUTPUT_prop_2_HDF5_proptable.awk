# from a cleaned list of variables defining the GALAXY_OUPUT C struct, create a .h file containing the
# LGalaxy and (if relevant) MOMAF output structs.
BEGIN{
}
{
	line=$0
	split(line,fields)
	type=fields[1]
	if(type=="//")
		next

	fieldname=fields[2]
	if(type == "long" && fieldname == "long") 
    	fieldname=fields[3]
	ia=index(fieldname,";")
	fieldname=substr(fieldname,0,ia-1)
	ib=index(fieldname,"[")
	if(ib!=0)
		fieldname=substr(fieldname,0,ib-1)

    unitidx=match(line,"//")
    if(unitidx!=0){
		slice_line=substr(line,unitidx+2)
		descidx=match(slice_line,"//")
		unit=substr(slice_line,0,descidx-1)
		desc=substr(slice_line,descidx+2)
		print fieldname "," unit "," desc
	}
	else
		print fieldname ", ," 
    
}
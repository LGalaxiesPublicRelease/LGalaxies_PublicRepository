# from a cleaned list of variables defining the GALAXY_OUPUT C struct, create a .h file containing the
# LGalaxy and (if relevant) MOMAF output structs.
BEGIN{
print "struct LGalaxy {"
}
{
	print $0 ";"
}
END{
print "};"
}
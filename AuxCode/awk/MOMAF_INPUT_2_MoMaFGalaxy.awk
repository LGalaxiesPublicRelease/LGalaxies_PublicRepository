# from a cleaned list of variables defining the MOMAF_INPUT C struct, create a cleaned MOMAF output structs.
BEGIN{
print "struct MoMaFGalaxy {"
}
{
	print $0 ";"
}
END{
print "};"
}
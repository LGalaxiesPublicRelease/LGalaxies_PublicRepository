/*------------------------------------------------------------------------------
  c interface to fortran subroutines in findtop.s.f
*/

void findtop(), findbot(), findtpa(), findbta(),
     finitop(), finibot(), finitpa(), finibta();

void findtop_(), findbot_(), findtpa_(), findbta_(),
     finitop_(), finibot_(), finitpa_(), finibta_();

void findtop(a, na, iord, nb)
int na, nb;
double a[];
int iord[];
{
    findtop_(a, &na, iord, &nb);
}

void findbot(a, na, iord, nb)
int na, nb;
double a[];
int iord[];
{
    findbot_(a, &na, iord, &nb);
}

void findtpa(a, na, iord, nb)
int na, nb;
double a[];
int iord[];
{
    findtpa_(a, &na, iord, &nb);
}

void findbta(a, na, iord, nb)
int na, nb;
double a[];
int iord[];
{
    findbta_(a, &na, iord, &nb);
}

void finitop(a, na, iord, nb)
int na, nb;
double a[];
int iord[];
{
    finitop_(a, &na, iord, &nb);
}

void finibot(a, na, iord, nb)
int na, nb;
double a[];
int iord[];
{
    finibot_(a, &na, iord, &nb);
}

void finitpa(a, na, iord, nb)
int na, nb;
double a[];
int iord[];
{
    finitpa_(a, &na, iord, &nb);
}

void finibta(a, na, iord, nb)
int na, nb;
double a[];
int iord[];
{
    finibta_(a, &na, iord, &nb);
}

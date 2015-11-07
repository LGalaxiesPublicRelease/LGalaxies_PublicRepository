#include <stdio.h>

struct metals
{
  float type1a;
  float type2;
  float agb;
};

struct metals metals_add(struct metals m1, 
	       struct metals m2,
	       float fraction);
struct metals metals_init();
void metals_print(char s[],struct metals m);
float metals_total(struct metals m);

int main()
{
  struct metals m1,m2;
  m1=metals_init();
  m2=metals_init();
  m1.type1a=1.;
  m2.agb=2.;
  metals_print("m1",m1);
  metals_print("m2",m2);
  m1=metals_add(m1,m2,0.1);
  metals_print("m1",m1);
  metals_print("m2",m2);

  return(0);
}

struct metals metals_add(struct metals m1, 
			 struct metals m2,
			 float fraction)
{
  struct metals m;
  m.type1a=m1.type1a+fraction*m2.type1a;
  m.type2=m1.type2+fraction*m2.type2;
  m.agb=m1.agb+fraction*m2.agb;

  return(m);
}

struct metals metals_init() 
{
  struct metals m;
  m.type1a=0.;
  m.type2=0.;
  m.agb=0.;
  return(m);
}

void metals_print(char s[],struct metals m) 
{
  printf("%s.type1a=%f\n",s,m.type1a);
  printf("%s.type2 =%f\n",s,m.type2);
  printf("%s.agb   =%f\n",s,m.agb);
  return;
}

float metals_total(struct metals m);
{
  return(m.type1a+m.type1b+m.agb);
}

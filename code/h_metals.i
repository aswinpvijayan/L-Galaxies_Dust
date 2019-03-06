# 1 "./code/h_metals.h"
# 1 "/lustre/scratch/astro/ap629/sc558/Dust_Rob_copy//"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "./code/h_metals.h"

struct metals
{
  float type1a;
  float type2;
  float agb;
};


//Individual element histories:
struct elements
{
  float H;
  float He;

  float Cb; //NOTE: Carbon (C) is stored as Cb here
  float N;

  float O;

  float Ne;

  float Mg;

  float Si;
  float S;
  float Ca;

  float Fe;
};
//Number of chemical elements tracked:
# 44 "./code/h_metals.h"
struct DustRates
{
  float AGB;
  float SNII;
  float SNIA;
  float GROW;
  float DEST;
};

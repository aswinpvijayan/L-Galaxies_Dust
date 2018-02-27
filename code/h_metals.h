#ifdef DETAILED_METALS_AND_MASS_RETURN
struct metals
{
  float type1a;
  float type2;
  float agb;
};

#ifdef INDIVIDUAL_ELEMENTS
//Individual element histories:
struct elements
{
  float H;
  float He;
#ifndef MAINELEMENTS
  float Cb; //NOTE: Carbon (C) is stored as Cb here
  float N;
#endif
  float O;
#ifndef MAINELEMENTS
  float Ne;
#endif
  float Mg;
#ifndef MAINELEMENTS
  float Si;
  float S;
  float Ca;
#endif
  float Fe;
};
//Number of chemical elements tracked:
#ifndef MAINELEMENTS
#define NUM_ELEMENTS 11 //All: [H][He][C][N][O][Ne][Mg][Si][S][Ca][Fe]
#else
#define NUM_ELEMENTS 5 //Only [H][He][O][Mg][Fe]
#endif

#endif //INDIVIDUAL_ELEMENTS
#endif //DETAILED_METALS_AND_MASS_RETURN

#ifdef DETAILED_DUST 			//Scott 15/11/2015
#define NDUST 5 // Needed for HDF5 table creation - should match struct below
struct DustRates
{
	float AGB;
	float SNII;
	float SNIA;
	float GROW;
	float DEST;
};
#endif //DETAILED_DUST

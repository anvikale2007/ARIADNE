#define NUM_DRUGS 8
#define MAX_STEP 12
const char* drugNames[] = {"No drug","Carbamazepine","Lamotrigine","Valproate","Cenobamate","Pregabalin","Levetiracetam","Ethosuximide"};
const uint8_t recruit[NUM_DRUGS][10] = {
  {3,255,10,6,9,5,10,6,7,6},
  {3,255,255,7,255,6,255,8,9,8},
  {3,255,255,7,255,6,255,7,9,8},
  {3,255,255,7,255,6,255,8,10,8},
  {3,255,255,7,255,6,255,8,10,9},
  {2,255,255,7,255,6,255,8,9,8},
  {2,255,255,6,255,6,255,7,9,8},
  {3,255,10,6,9,5,10,6,7,6},
};

#ifndef NewPulseShapes_h
#define NewPulseShapes_h 1

#include <typeinfo>
#include <vector>

class NewPulseShapes {
 public:
  NewPulseShapes();
  ~NewPulseShapes();

  void init();
  void configurePulseShapes();

  float getPulseFrac(float fC, float time, int ts) const;
  float getPulseFracNorm(float fC, float time) const;
  //float getTimeSlew(float fC, int ts) const;

 private:
  float loThresh;
  float hiThresh;
  float flip[10];
  float par0[10][2];
  float par1[10][2];
  float par2[10][2];
  float par3[10][2];

  float tpar0[10];
  float tpar1[10];
  float tpar2[10];
};

#endif

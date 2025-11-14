#include <Arduino.h>
#include "LD6002.h"
LD6002 radar(Serial1);

void setup()
{
  Serial.begin(115200);
  Serial1.begin(1382400, SERIAL_8N1, 21, 20);

}

float lastHeartRate = 0;
float lastBreathRate = 0;
float lastDistance = 0;
void loop()
{
  radar.update();

  if (radar.hasNewHeartRate())
  {
   float heartRateMain = radar.getHeartRate();
   if ((heartRateMain != lastHeartRate) && (heartRateMain > 0))
   {
    Serial.printf("Heart Rate: %.2f\n", heartRateMain);
   }
   lastHeartRate = heartRateMain;
   radar.clearHeartRateFlag();
  }

  if (radar.hasNewBreathRate())
  {
   float breathRateMain = radar.getBreathRate();
   if ((breathRateMain != lastBreathRate) && (breathRateMain > 0))
   {
    Serial.printf("Breath Rate: %.2f\n", breathRateMain);
   }
   lastBreathRate = breathRateMain;
   radar.clearBreathRateFlag();
  }

  if (radar.hasNewDistance())
  {
   float distanceMain = radar.getDistance();
   if ((distanceMain != lastDistance) && (distanceMain > 0))
   {
    Serial.printf("Distance: %.2f\n", distanceMain);
   }
   lastDistance = distanceMain;
   radar.clearDistanceFlag();
  }
}

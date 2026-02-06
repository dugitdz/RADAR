#include <Arduino.h>
#include <stdint.h>

void setup() {
  Serial.begin(1382400);
  Serial1.begin(1382400, SERIAL_8N1, 21, 20);
  Serial1.setTimeout(60);
}

void loop() {
  while (Serial1.available()) {
    uint8_t b = (uint8_t)Serial1.read();
    Serial.write(b);   // manda o byte CRU pro PC
  }
}

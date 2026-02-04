#include <Arduino.h>
#include "LD6002.h"
#include <stdint.h>

// ---------- CÓDIGO 1: LD6002 (x15, x14) ----------

LD6002 radar(Serial1);

float lastHeartRate  = 0;
float lastBreathRate = 0;
float lastDistance   = 0;

// ---------- CÓDIGO 2: LEITURA DIRETA (x13) ----------

float leitor_bytes_phase(uint8_t *aux);

// variáveis do parser x13 (mantidas como no seu código)
uint8_t prox_byte = 0;
static uint8_t buffer[32];
int len = 0;
float heart_phase, total_phase, breath_phase, dist;
static unsigned long t_phase = 0, t_dist = 0;
static bool b_phase = false, b_dist = false;
static const unsigned long MAX_dt = 100;

// criador de CSV: exatamente o mesmo print que você já fazia
void escreveCSV(unsigned long t_out,
                float total_phase,
                float breath_phase,
                float heart_phase,
                float dist)
{
  Serial.print(t_out);       Serial.print(",");
  Serial.print(total_phase); Serial.print(",");
  Serial.print(breath_phase);Serial.print(",");
  Serial.print(heart_phase); Serial.print(",");
  Serial.println(dist);
}

void setup()
{
  // Mantendo o baud do SEGUNDO código (você pode mudar se quiser)
  Serial.begin(1382400);

  // Serial1 igual em ambos os códigos
  Serial1.begin(1382400, SERIAL_8N1, 21, 20);
  Serial1.setTimeout(60);

  // Cabeçalho CSV do segundo código
  Serial.println("tempo,total_phase,breath_phase,heart_phase,dist");
}

void loop()
{
  // ---------- PARTE 1: código do LD6002 (x15, x14) ----------

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

  // ---------- PARTE 2: código do parser x13 (mantido) ----------

  while (Serial1.available())
  {
    uint8_t s = Serial1.read();
    if (s == 0x0A)
    {
      if (Serial1.readBytes(&prox_byte, 1) != 1) continue;

      switch (prox_byte)
      {
      case 0x13:
        len = 14;
        if (len > (int)sizeof(buffer)) break;
        buffer[0] = prox_byte;
        {
          size_t tam = Serial1.readBytes(&buffer[1], len - 1);
          if (tam != len - 1) break;
        }
        t_phase = millis();
        total_phase  = leitor_bytes_phase(&buffer[2]);
        breath_phase = leitor_bytes_phase(&buffer[6]);
        heart_phase  = leitor_bytes_phase(&buffer[10]);
        b_phase = true;
        if (b_dist)
        {
          unsigned long dt = (t_phase > t_dist) ? (t_phase - t_dist) : (t_dist - t_phase);
          if (dt <= MAX_dt)
          {
            unsigned long t_out = (t_phase > t_dist) ? t_phase : t_dist;
            // AQUI usamos o criador de CSV (mesmo print que você já tinha)
            escreveCSV(t_out, total_phase, breath_phase, heart_phase, dist);
            b_dist = false;
            b_phase = false;
          }
        }
        break;

      case 0x16:
        len = 10;
        if (len > (int)sizeof(buffer)) break;
        buffer[0] = prox_byte;
        {
          size_t tam = Serial1.readBytes(&buffer[1], len - 1);
          if (tam != len - 1) break;
        }
        t_dist = millis();
        dist = leitor_bytes_phase(&buffer[6]);
        b_dist = true;
        if (b_phase)
        {
          unsigned long dt = (t_phase > t_dist) ? (t_phase - t_dist) : (t_dist - t_phase);
          if (dt <= MAX_dt)
          {
            unsigned long t_out = (t_phase > t_dist) ? t_phase : t_dist;
            // Mesmo CSV
            escreveCSV(t_out, total_phase, breath_phase, heart_phase, dist);
            b_dist = false;
            b_phase = false;
          }
        }
        break;

      case 0x14:
      default:
        break;
      }
    }
  }
}

// função de conversão de bytes -> float (mantida igual)
float leitor_bytes_phase(uint8_t *aux)
{
  uint32_t u = (uint32_t)aux[0]
             | ((uint32_t)aux[1] << 8)
             | ((uint32_t)aux[2] << 16)
             | ((uint32_t)aux[3] << 24);
  float f;
  memcpy(&f, &u, 4);
  return f;
}

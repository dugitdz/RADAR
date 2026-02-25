#include <Arduino.h>
#include <stdint.h>
#include <string.h>

#include <BLEDevice.h>
#include <BLEServer.h>
#include <BLEUtils.h>
#include <BLE2902.h>

static const uint32_t RADAR_BAUD = 1382400;
static const int RADAR_RX = 21;
static const int RADAR_TX = 20;

static const char *BLE_NAME     = "ESP32C3_IC";
static const char *SERVICE_UUID = "ABCD";
static const char *CHAR_UUID    = "1234";

BLECharacteristic *g_char = nullptr;
static volatile bool g_connected = false;

class MyServerCallbacks : public BLEServerCallbacks {
  void onConnect(BLEServer *pServer) override { g_connected = true; }
  void onDisconnect(BLEServer *pServer) override {
    g_connected = false;
    BLEDevice::startAdvertising();
  }
};

static inline void pack_float_le(float x, uint8_t out4[4]) {
  memcpy(out4, &x, 4);
}

void setup() {

  // Radar
  Serial1.begin(RADAR_BAUD, SERIAL_8N1, RADAR_RX, RADAR_TX);
  Serial.begin(115200);
  Serial1.setTimeout(60);

  // BLE
  BLEDevice::init(BLE_NAME);
  BLEServer *server = BLEDevice::createServer();
  server->setCallbacks(new MyServerCallbacks());

  BLEService *service = server->createService(SERVICE_UUID);

  g_char = service->createCharacteristic(
    CHAR_UUID,
    BLECharacteristic::PROPERTY_READ |
    BLECharacteristic::PROPERTY_WRITE |
    BLECharacteristic::PROPERTY_NOTIFY
  );
  g_char->addDescriptor(new BLE2902());

  service->start();

  BLEAdvertising *adv = BLEDevice::getAdvertising();
  adv->addServiceUUID(SERVICE_UUID);
  adv->start();
}

void loop() {

  static const int REST_A13 = 14;
  static const int REST_A14 = 6;
  static const int REST_A15 = 6;
  static const int REST_A16 = 10;

  static uint8_t restbuf[32];

  static uint8_t out_pkt[15 + 4];

  while (Serial1.available()) {
    uint8_t b = (uint8_t)Serial1.read();
    if (b != 0x0A) continue; // 1o byte do TYPE (0x0A)

    uint8_t type2;
    if (Serial1.readBytes(&type2, 1) != 1) continue; // 2o byte do TYPE (0x13/14/15/16)

    int need = 0;
    switch (type2) {
      case 0x13: need = REST_A13; break;
      case 0x14: need = REST_A14; break;
      case 0x15: need = REST_A15; break;
      case 0x16: need = REST_A16; break;
      default:
        continue;
    }

    if (need > (int)sizeof(restbuf)) continue;

    size_t got = Serial1.readBytes(restbuf, (size_t)need);
    if (got != (size_t)need) continue;

    if (type2 == 0x13) {

      out_pkt[0] = 0x13;
      memcpy(&out_pkt[1], restbuf, REST_A13);

      float t_ms = (float)millis();
      Serial.print("\nTimestamp (ms): ");
      Serial.println(t_ms);
      uint8_t t4[4];
      pack_float_le(t_ms, t4);

      memcpy(&out_pkt[15], t4, 4);

      // BLE notify
      if (g_connected && g_char) {
        g_char->setValue(out_pkt, sizeof(out_pkt));
        g_char->notify();
      }
    }
  }
}

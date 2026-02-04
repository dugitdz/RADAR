import asyncio
import struct
from bleak import BleakScanner, BleakClient

DEVICE_NAME = "ESP32C3_IC"
CHAR_UUID   = "00001234-0000-1000-8000-00805f9b34fb"

def read_float_le(b4: bytes) -> float:
    # equivalente ao teu leitor_bytes_phase (little-endian)
    u = (b4[0]
         | (b4[1] << 8)
         | (b4[2] << 16)
         | (b4[3] << 24)) & 0xFFFFFFFF
    return struct.unpack("<f", u.to_bytes(4, "little"))[0]

def on_notify(sender: int, data: bytearray):
    # Esperado: 19 bytes = 15 (A13 bruto) + 4 (timestamp float)
    if len(data) != 19:
        print(f"[WARN] len={len(data)} hex={data.hex()}")
        return

    if data[0] != 0x13:
        print(f"[WARN] primeiro byte != 0x13: {data[0]:02X} | hex={data.hex()}")
        return

    total_phase  = read_float_le(data[2:6])
    breath_phase = read_float_le(data[6:10])
    heart_phase  = read_float_le(data[10:14])
    t_ms         = read_float_le(data[15:19])

    print(f"t_ms={t_ms:9.2f} | total={total_phase: .6f} | breath={breath_phase: .6f} | heart={heart_phase: .6f}")

async def main():
    print("Scanning BLE...")
    devices = await BleakScanner.discover(timeout=5.0)

    target = None
    for d in devices:
        if d.name == DEVICE_NAME:
            target = d
            break

    if not target:
        print(f"Não achei '{DEVICE_NAME}'. Vi:")
        for d in devices:
            print(" ", d.address, d.name)
        return

    print("Conectando em:", target.address, target.name)

    async with BleakClient(target.address) as client:
        if not client.is_connected:
            print("Falhou conectar.")
            return

        print("Conectado. start_notify...")
        await client.start_notify(CHAR_UUID, on_notify)

        print("Recebendo (Ctrl+C pra sair)...")
        try:
            while True:
                await asyncio.sleep(1.0)
        except KeyboardInterrupt:
            pass

        await client.stop_notify(CHAR_UUID)
        print("Parou.")

if __name__ == "__main__":
    asyncio.run(main())

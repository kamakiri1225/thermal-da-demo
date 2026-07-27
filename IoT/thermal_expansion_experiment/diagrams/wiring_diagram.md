# 配線図（結線表は docs/04_electrical_wiring.md 参照）

## Mermaidによる結線イメージ

```mermaid
flowchart TB
    ESP32["ESP32-DevKitC-32E"]

    ESP32 -- "GPIO18 SCK" --> BUS_SCK["SPIバス: SCK"]
    ESP32 -- "GPIO19 MISO" --> BUS_MISO["SPIバス: MISO"]
    ESP32 -- "GPIO23 MOSI" --> BUS_MOSI["SPIバス: MOSI"]

    BUS_SCK --> M1["MAX31856 #1"]
    BUS_SCK --> M2["MAX31856 #2"]
    BUS_SCK --> M3["MAX31856 #3"]
    BUS_MISO --> M1
    BUS_MISO --> M2
    BUS_MISO --> M3
    BUS_MOSI --> M1
    BUS_MOSI --> M2
    BUS_MOSI --> M3

    ESP32 -- "GPIO25 CS1" --> M1
    ESP32 -- "GPIO26 CS2" --> M2
    ESP32 -- "GPIO27 CS3" --> M3

    M1 --- TC1["K型熱電対 T1"]
    M2 --- TC2["K型熱電対 T2"]
    M3 --- TC3["K型熱電対 T3"]

    ESP32 -- "GPIO34 DATA(入力専用)" --> IND["デジマチックインジケータ"]
    ESP32 -- "GPIO35 CLOCK(入力専用)" --> IND
    ESP32 -- "GND共通" --- IND

    ESP32 -- "3V3" --> PWR["電源レール 3.3V"]
    PWR --> M1
    PWR --> M2
    PWR --> M3
```

## ASCII結線図（詳細はdocs/04_electrical_wiring.md）

```
                        ┌──────────────┐
   K型熱電対T1 ──T+/T−→│  MAX31856 #1 │
   (加熱側)             │  CS=GPIO25   │──┐
                        └──────────────┘  │
                        ┌──────────────┐  │   SPIバス共通
   K型熱電対T2 ──T+/T−→│  MAX31856 #2 │──┤   SCK  = GPIO18
   (中央)                │  CS=GPIO26   │  │   MISO = GPIO19
                        └──────────────┘  │   MOSI = GPIO23
                        ┌──────────────┐  │
   K型熱電対T3 ──T+/T−→│  MAX31856 #3 │──┘
   (自由端側)            │  CS=GPIO27   │
                        └──────────────┘

   デジマチックインジケータ (電池駆動、ESP32から給電しない)
        DATA  ───────────→ GPIO34 (入力専用)
        CLOCK ───────────→ GPIO35 (入力専用)
        GND   ───────────── GND共通

   すべてのMAX31856: VIN→3V3, GND→GND
```

## 注意（再掲）

- Digimatic側の配線・保護回路は`docs/05_digimatic_protocol.md`の安全確認手順を先に実施してから通電すること。
- AC100V配線はこの図には含めない（別系統、`docs/04_electrical_wiring.md`4.5節・`docs/10_safety.md`参照）。

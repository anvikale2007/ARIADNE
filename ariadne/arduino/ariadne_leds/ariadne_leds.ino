// ARIADNE mini-brain demo: WS2812B (NeoPixel) strip, 12 LEDs.
// Receives a line like  S:001221000110\n  over serial (115200 baud) from led.py.
// 0 = green (quiet), 1 = amber (partly recruited), 2 = red (seizing).
#include <Adafruit_NeoPixel.h>
#define PIN 6
#define N_LEDS 12
Adafruit_NeoPixel strip(N_LEDS, PIN, NEO_GRB + NEO_KHZ800);

void setup() {
  Serial.begin(115200);
  strip.begin(); strip.setBrightness(80);
  for (int i = 0; i < N_LEDS; i++) strip.setPixelColor(i, strip.Color(0, 60, 0));
  strip.show();
}

void loop() {
  if (Serial.available()) {
    String line = Serial.readStringUntil('\n');
    if (line.startsWith("S:")) {
      for (int i = 0; i < N_LEDS && i + 2 < (int)line.length(); i++) {
        char c = line.charAt(i + 2);
        uint32_t col = strip.Color(0, 80, 0);                    // green
        if (c == '1') col = strip.Color(120, 70, 0);             // amber
        if (c == '2') col = strip.Color(150, 0, 0);              // red
        strip.setPixelColor(i, col);
      }
      strip.show();
    }
  }
}

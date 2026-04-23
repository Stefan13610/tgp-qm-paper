# Q7: Dekoherencja z dynamiki pola TGP

## Problem

Wyjasnic przejscie kwantowo-klasyczne (dekoherencje) jako naturalna
konsekwencje dynamiki pola TGP, bez postulowania oddzielnego mechanizmu.

## Mechanizm -- TRZY DROGI DO DEKOHERENCJI

### Droga 1: Tlumienie hbar(Phi)
- hbar(Phi) = hbar_0 * sqrt(Phi_0/Phi)
- Gestzy substrat (Phi >> Phi_0) => hbar -> 0 => klasyczny
- Tempo dekoherencji: Gamma_1 ~ eps/(1+eps)

### Droga 2: Nieliniowe mieszanie modow (z Q3)
- Superpozycja liniowa dla eps << 1
- Korekta NL: delta_A/A ~ 0.392*eps
- Mieszanie modow (k=1 -> k=2) niszczy koherencje fazowa
- Tempo: Gamma_2 ~ eps^2

### Droga 3: Srodowiskowa back-reakcja (z Q4)
- N_env solitonow w odleglosci D z amplituda A_env
- Kazdy rozprzsza faze: delta_phi ~ chi*A_env/D
- Bladzenie losowe: <delta_phi^2> ~ N_env * chi^2 * A_env^2/D^2
- Tempo: Gamma_3 ~ eps * (A/D)

WSZYSTKIE TRZY sa EMERGENTNE z dynamiki TGP, nie postulowane.

## Wyniki (2026-04-15) -- q7_decoherence.py (8/8 PASS)

### hbar(Phi) i limit klasyczny
- hbar(Phi) = hbar_0 * sqrt(Phi_0/Phi): zweryfikowany
- hbar(1e10) = 1e-5 (klasyczny), hbar(0.01) = 10 (kwantowy)
- Przy Ziemi: delta_hbar/hbar ~ -3.5e-10

### Dlugosc koherencji
- L_coh = 2*Phi / |grad(Phi)| = 2*r^2/A (blisko solitonu)
- Daleko od solitonu: L_coh >> r (koherentny)
- Blisko solitonu: L_coh ~ r (dekoherentny)

### Tempo dekoherencji NL
- Gamma_NL ~ N^2 dokladnie (skalowanie potwierdzone)
- N_class = D/(0.39*A) ~ 2551 dla typowego ukladu
- Zgadza sie z Q3: limit klasyczny przy N ~ D/A_tail

### Dekoherencja srodowiskowa
- Powietrze (STP): Gamma = 2e+28 (natychmiastowa dekoherencja)
- Pustka kosmiczna: Gamma = 8e-15 (koherencja zachowana)
- Stosunek: 2.4e+42 (!)

### Macierz gestosci
- Elementy pozadiagonalne: rho_ij(t) = rho_ij(0) * exp(-Gamma*t)
- Tr(rho) = 1 zachowane (populacje nienaruszone)
- Czystosc: 1.000 (czysty) -> 0.500 (mieszany)

### Trzy drogi -- ujednolicenie
- Droga 1 dominuje (Gamma ~ eps), drogi 2,3 ~ eps^2
- Wszystkie zanikaja gdy eps -> 0 (rezsim kwantowy odzyskany)
- Hierarchia: hbar-tlumienie >> NL-mieszanie >> back-reakcja

### Darwinizm kwantowy
- Informacja wzajemna saturuje po 23.3% srodowiska
- Redundancja R = 4.3x (wielokrotne kopie informacji)
- TGP wyjasnienie: ogon solitonu odciska faze na CALYM substracie

## Lancuch derywacji

```
Dynamika pola TGP
  |
  +---> Droga 1: hbar(Phi) = hbar_0*sqrt(Phi_0/Phi)
  |     Gesty substrat => hbar -> 0 => klasyczny
  |
  +---> Droga 2: NL mieszanie modow (Q3)
  |     eps duze => superpozycja lamie sie => klasyczny
  |
  +---> Droga 3: Srodowiskowa back-reakcja (Q4)
  |     N_env solitonow rozpraszaja faze => dekoherencja
  |
  v
DEKOHERENCJA jest EMERGENTNA, nie postulowana
Darwinizm kwantowy z geometrycznego nakladania ogonow
```

## Predykcje testowalne

1. Zmiennosc hbar(Phi): eksperymenty kwantowe na roznych wysokosciach
   delta_hbar/hbar ~ -GM/(2*rc^2) ~ -3.5e-10 (Ziemia)
2. Tempo dekoherencji TGP vs standardowa QM (rozne skalowanie z N)
3. Dlugosc koherencji: L_coh ~ r^2/A (przestrzenna zmiennosc)
4. Modyfikacja BEC: T_BEC ~ hbar(Phi)^2

## STATUS: Q7 ZAMKNIETE

- [x] hbar(Phi) forma funkcyjna i limit klasyczny
- [x] Dlugosc koherencji z gradientu hbar
- [x] NL tempo dekoherencji ~ eps^2
- [x] Srodowiskowa dekoherencja z back-reakcji
- [x] Ewolucja macierzy gestosci
- [x] Trzy drogi ujednolicone
- [x] Darwinizm kwantowy z TGP

## Pliki

| Plik | Opis | Status |
|------|------|--------|
| q7_decoherence.py | hbar(Phi), koherencja, NL, srodowisko, darwinizm | 8/8 PASS |

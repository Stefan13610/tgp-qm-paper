# Q3: Superpozycja z liniowosci slabego pola

## Problem

Wyjasnic superpozycje kwantowa jako konsekwencje rownan TGP,
a nie jako fundamentalny postulat.

## Mechanizm

Pelne ODE substratu: g'' + (1/g)g'^2 + (2/r)g' = 1 - g

W poblizu vacuum (g ~ 1), niech g = 1 + eps*f:
  f'' + (2/r)f' + f = -eps*f'^2/(1+eps*f)

Dla eps -> 0: **f'' + (2/r)f' + f = 0** -- rownanie LINIOWE!
Rozwiazania sie superpozycjonuja: kazda kombinacja liniowa jest rozwiazaniem.

Nieliniowa korekta RHS = -eps*f'^2/(1+eps*f) lamie superpozycje
dla silnych pol (duze eps).

## Wyniki (2026-04-15) -- q3_superposition.py (7/7 PASS)

### Walidacja linearyzacji
- eps=0.01: blad < 0.4% (quasi-ideal)
- eps=0.13 (elektron): korekta ~5% (quasi-liniowy rezim)
- eps=0.5: korekta ~20% (silna nieliniownosc)

### Test superpozycji
- Blad superpozycji ~ eps^0.98 (R^2=0.999) -- dokladnie liniowy
- Dla eps < 0.01: blad < 0.1%
- Superpozycja jest EMERGENTNA, nie fundamentalna

### Korekta nieliniowa
- dA/A = 0.392 * eps^1.01 (R^2=0.9999)
- eps_crit(10%) = 0.26
- eps_crit(50%) = 1.27
- Skalowanie dokladnie liniowe w eps

### Druga harmoniczna (sygnatura nieliniowosci)
- A_2/A_1 = 0.03% dla elektronu (g0=0.869)
- A_2/A_1 rosnie z amplituda solitonu
- Predykcja testowalna: drugie harmoniczne w rozpraszaniu

### Limit klasyczny z TGP
- N solitonow o amplitudzie A_tail w odleglosci D:
  eps_eff ~ N*A_tail/D
- Superpozycja lamie sie gdy eps_eff >> 1
- N_class(atom) ~ 800 czastek
- Obiekty makroskopowe (N >> 10^23): calkowicie klasyczne

### Dekoherencja
- Nieliniowe mieszanie modow (k=1 -> k=2) rozprasza informacje fazowa
- Tempo dekoherencji ~ eps^2 ~ N^2*A^2/D^2
- Srodowiskowa dekoherencja WYNIKA z dynamiki TGP

## Lancuch derywacji

```
TGP ODE -> linearyzacja przy vacuum -> rownanie LINIOWE
-> superpozycja rozwiazan -> SUPERPOZYCJA KWANTOWA
-> korekty nieliniowe przy silnym polu -> LIMIT KLASYCZNY
-> mieszanie modow z NL -> DEKOHERENCJA
```

## STATUS: Q3 ZAMKNIETE ✓

- [x] Walidacja linearyzacji (blad ~ eps)
- [x] Test superpozycji przeszedl dla slabych pol
- [x] Skalowanie korekcji nieliniowej ustalone
- [x] Limit klasyczny z N-solitonowej nieliniowosci
- [x] Polaczenie z dekoherencja zidentyfikowane

## Pliki

| Plik | Opis | Status |
|------|------|--------|
| q3_superposition.py | Linearyzacja, superpozycja, NL korekty | 7/7 PASS |

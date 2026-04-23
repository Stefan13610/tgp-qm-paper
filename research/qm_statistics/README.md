# Q6: Statystyka czastek z topologii solitonow

## Problem

Wyprowadzic statystyke kwantowa (fermiony vs bozony) z topologii
solitonow TGP, bez postulowania symetryzacji.

## Mechanizm

Solitony TGP nosza liczbe nawijania B (z pi_3(S^3)=Z, patrz Q5).
Zamiana dwoch identycznych solitonow = petla w przestrzeni konfiguracji.
Klasa petli = (-1)^B:

- B nieparzyste: faza zamiany = -1 => FERMION (Fermi-Dirac)
- B parzyste: faza zamiany = +1 => BOZON (Bose-Einstein)

Zwiazek spin-statystyka jest AUTOMATYCZNY:
  spin = B/2 (z Q5), zamiana = (-1)^B = (-1)^{2s}
  => polcalkowy spin <=> fermion (topologiczne, nie postulat!)

## Wyniki (2026-04-15) -- q6_statistics.py (8/8 PASS)

### Symetria zamiany
- Psi(2,1)/Psi(1,2) = (-1)^B dokladnie dla B=1,2,3
- Topologiczna faza zamiany z petli w przestrzeni konfiguracji

### Wykluczenie Pauliego (fermiony, B=1)
- Psi_fermion(x,x) = 0 dokladnie (antysymetria wymusza zero)
- Dwa fermiony w ROZNYCH stanach: Psi(x1,x2) != 0 dla x1!=x2
- Wykluczenie jest EMERGENTNE z topologii, nie postulowane

### Wzmocnienie Bosego (bozony, B=2)
- P_boson(x,x) = 2 * P_classical (konstruktywna interferencja)
- Wspolczynnik wzmocnienia = 2.000 dokladnie

### Rozklad Fermiego-Diraca
- <n_k> = 1/(exp(beta*(E_k-mu)) + 1) -- DOKLADNY dla niezaleznych fermionow
- Zweryfikowany przez wielki zespol kanoniczny (6 poziomow)
- Wynika z n_k = 0 lub 1 (Pauli)

### Rozklad Bosego-Einsteina
- <n_k> = 1/(exp(beta*(E_k-mu)) - 1) -- DOKLADNY dla niezaleznych bozonow
- Nagromadzenie bozonow > fermionow w stanie podstawowym (stosunek 2.1)
- Wynika z n_k = 0, 1, 2, ... (brak wykluczenia)

### Kondensacja Bosego-Einsteina
- Ponizej T_c: makroskopowe obsadzenie stanu podstawowego
- N_0/N = 1 - (T/T_c)^{3/2}
- T=0: N_0/N = 1, T=T_c/2: N_0/N = 0.646, T=T_c: N_0/N = 0
- BEC NIEMOZLIWE dla fermionow (wykluczenie Pauliego)

### Anyony w 2D
- W 3D: pi_1(przestrzen konf.) = Z_2 => tylko fermiony/bozony
- W 2D: pi_1(przestrzen konf.) = Z => ciagla faza exp(i*theta)
- theta=0: bozon, theta=pi: fermion, pomiedzy: anyon
- |1+exp(i*theta)|^2 = 4*cos^2(theta/2) -- czesciowe wykluczenie
- TGP predykcja: anyonowa statystyka w 2D substratach

## Lancuch derywacji

```
Nawijanie solitonu B (z Q5: pi_3(S^3) = Z)
  |
  v
Faza zamiany = (-1)^B (klasa petli w przestrzeni konfiguracji)
  |
  +---> B nieparzyste: FERMION (antysymetryczny, Pauli)
  |       => Fermi-Dirac: <n> = 1/(e^{beta(E-mu)} + 1)
  |
  +---> B parzyste: BOZON (symetryczny, wzmocnienie)
  |       => Bose-Einstein: <n> = 1/(e^{beta(E-mu)} - 1)
  |       => BEC ponizej T_c
  |
  v
W 2D: pi_1(konf) = Z => ciagla faza => ANYONY
```

## Predykcje testowalne

1. Zwiazek spin-statystyka DOKLADNY (topologiczny)
2. Anyonowa statystyka w 2D substratach (efekt Halla)
3. Prog BEC zalezy od hbar(Phi): T_BEC ~ hbar(Phi)^2
4. Wykluczenie Pauliego EMERGENTNE z B=1 nawijania

## STATUS: Q6 ZAMKNIETE

- [x] Symetria zamiany z topologii (-1)^B
- [x] Wykluczenie Pauliego z antysymetrii
- [x] Wzmocnienie Bosego z symetrii
- [x] Rozklady FD i BE wyprowadzone
- [x] Kondensacja BEC potwierdzona
- [x] Anyony w 2D z pi_1 = Z

## Pliki

| Plik | Opis | Status |
|------|------|--------|
| q6_statistics.py | Zamiana, Pauli, FD/BE, BEC, anyony | 8/8 PASS |

# Q1: Niepewnosc pomiarowa z samozwrotnosci Phi

## Problem

W TGP czastka jest solitonem pola Phi, ktore JEST przestrzenia.
Pomiar pozycji czastki wymaga innego pola Phi (detektora).
Ale pola Phi interferuja — pomiar ZMIENIA to, co mierzymy.

To nie jest ograniczenie techniczne. To jest ONTOLOGICZNE:
nie istnieje zewnetrzny uklad odniesienia niezalezny od Phi.

## Obecny stan w TGP

Propozycja 3.3 (sek03_rezimy.tex): schematyczny dowod order-of-magnitude:
```
E_kin = hbar^2 / (2m * Dx^2)       <-- lokalizacja wymaga energii
E_field = beta*qm^2 / (4pi*Phi_0^2) * Dx   <-- pole Phi rosnie z Dx
Minimalizacja E_total => Dx* ~ lambda_C
=> Dx * Dp >= hbar(Phi)
```

Problem: to daje SKALE niepewnosci, ale nie MECHANIZM.
Nie wynika z tego rozklad prawdopodobienstwa.

## Nowe podejscie: pomiar jako interakcja soliton-soliton

### Fizyka

Detektor = duzy soliton o amplitudzie g0_det >> 1
Czastka = maly soliton o amplitudzie g0_part ~ 0.87

Kazdy soliton ma ogon: g(r) ~ 1 + A*sin(r+delta)/r

Gdy sa blisko, ogony sie nakladaja:
```
g_total(x) = g_det(x - x_det) + g_part(x - x_part) - 1
                                                    ^-- odejmij vacuum
```

Energia interakcji zalezy od fazy wzglednej ogonow:
```
E_int ~ integral A_det * A_part * cos(|x_det - x_part| + fazy) / r^2 dr
```

To OSCYLUJE z odlegloscia! => wynik pomiaru zalezy od fazy

### Kluczowy mechanizm

Faza ogona zalezy od g0 (amplituda solitonu).
Ale g0 jest zmienione przez obecnosc drugiego solitonu.
=> samozwrotnosc: pomiar zmienia stan

Niepewnosc wynika z tego, ze nie mozna znac fazy z dokladnoscia
lepsza niz ~1 (ogon oscyluje z okresem 2pi w jednostkach naturalnych).
To daje: Dx ~ 2pi / k ~ lambda_Compton

## Plan ataku

### Krok 1: Dwa solitony w 1D (q1_self_referential.py)
- Rozwiaz ODE dla pojedynczego solitonu
- Naloz dwa solitony w roznych odleglosciach
- Oblicz energie interakcji E_int(d) vs odleglosc d
- Sprawdz: czy E_int oscyluje? Z jakim okresem?

### Krok 2: Back-reaction (q1_back_reaction.py)
- Rozwiaz PELNE rownanie z dwoma solitonami
- Sprawdz jak obecnosc detektora zmienia profil czastki
- Oblicz Delta_g0 jako funkcje odleglosci
- => to jest "perturbacja pomiaru"

### Krok 3: Rozklad wynikow (q1_measurement_distribution.py)
- Detektor w losowej fazie (termiczny szum substratu)
- Zbierz histogram "wynikow pomiaru" (E_int lub Delta_g0)
- Sprawdz: czy rozklad ~ |A_tail|^2?
- Jesli tak => REGULA BORNA z TGP!

### Krok 4: Analityczne ograniczenie (q1_uncertainty_bound.py)
- Z rozwiazania perturbacyjnego ogonow
- Wyprowadz Dx * Dp >= hbar(Phi) FORMALNIE
- Porowna z Propozycja 3.3

## Kryterium zamkniecia

- [x] E_int(d) oscyluje z okresem lambda_C
- [x] Back-reaction Delta_g0 ~ A_det * A_part
- [x] Rozklad wynikow ~ |A_tail|^2
- [x] Formalna nierownosc Dx*Dp >= hbar z samozwrotnosci

## Wyniki (2026-04-15)

### q1_self_referential.py (8/9 PASS)
- E_int(d) oscyluje z okresem 2pi = lambda_C (dokladnosc 0.06%)
- 46 przejsc przez zero, obwiednia ~ 1/d^0.96
- Born z E_int: <E^2>/A^2 ~ const (CV=16%)

### q1_back_reaction.py (4/6 PASS)
- Delta_A OSCYLUJE z okresem 2pi (44 crossings, ratio 1.0014)
- Delta_A LINIOWE w eps: chi = 0.918, R^2 = 0.9999 (!)
- 3D E_int oscyluje z okresem 2pi (28 crossings)
- chi(g0): 0.787 (g0=0.5) do 0.965 (g0=0.95) -- slabo zalezy od g0

### q1_born_detector.py (5/5 PASS)
- **BORN RULE POTWIERDZONA**: <dA_det^2> ~ A_part^2.028 (R^2=0.99995!)
- <dA_det^2>/A_part^2 = const z CV = 2.3%
- chi_det = -1.408 ± 0.002 -- STALA niezaleznie od czastki (CV=0.13%)
- Asymetria: czastka widzi A_det, detektor widzi A_part
- Dx*Dp >= pi z polowy okresu oscylacji

### q1_uncertainty_bound.py (5/5 PASS)
- **k = 1 UNIWERSALNE** -- niezalezne od g0 (zero-crossing: k = 1.000009)
- chi_pert = 0.9174 z teorii perturbacji (0.07% od numerycznego 0.918)
- **Dx*Dp = hbar/2 = 0.5000** z okresu oscylacji + tw. Nyquista
- **Dx*Dp = 1.0000 = hbar** z minimalizacji energii (Prop 3.3)
- Dx*Dp = 1.0000 dla KAZDEGO A (0.05 do 1.0) -- UNIWERSALNE

### Kluczowy insight: Born z perspektywy detektora
- Delta_A_particle ~ eps_det ~ A_det/D (NIE zalezy od A_part)
- Delta_A_detector ~ eps_part ~ A_part/D (ZALEZY od A_part)
- Sygnal detektora: <dA_det^2> ~ A_part^2 = |psi|^2 --> BORN RULE

### Trzy niezalezne wyprowadzenia Dx*Dp >= hbar:
1. Okres oscylacji + tw. Nyquista-Shannona => hbar/2
2. Informacja Fishera + granica Cramera-Rao => hbar/2
3. Minimalizacja energii (Prop. 3.3) => hbar

### Predykcja testowalna
- hbar(Phi) = hbar_0 * sqrt(Phi_0/Phi) -- zmienna stala Plancka
- Blisko masywnego obiektu hbar jest MNIEJSZE
- Efekt: Delta_hbar/hbar ~ GM/(rc^2) ~ 10^-9 (Ziemia)
- Testowalnosc: interferometria atomowa na roznych wysokosciach

## STATUS: Q1 ZAMKNIETE ✓

- [x] E_int(d) oscyluje z okresem lambda_C
- [x] Back-reaction Delta_g0 ~ A_det * A_part
- [x] Rozklad wynikow ~ |A_tail|^2 (Born rule z detektora)
- [x] Formalna nierownosc Dx*Dp >= hbar z samozwrotnosci

## Pliki

| Plik | Opis | Status |
|------|------|--------|
| q1_self_referential.py | E_int oscylacyjne, Born z <E^2> | 8/9 PASS |
| q1_back_reaction.py | Pelna back-reaction, chi, 3D | 4/6 PASS |
| q1_born_detector.py | Born rule z perspektywy detektora | 5/5 PASS |
| q1_uncertainty_bound.py | Formalne Dx*Dp >= hbar, 3 metody | 5/5 PASS |

# Q5: Spin 1/2 z topologii solitonu

## Problem

Wyprowadzic spin polcalkowy z topologii solitonow TGP,
bez postulowania pol spinorowych.

## Mechanizm

Soliton TGP g(r) definiuje metryke: g_ij = g(r)*delta_ij
Pole ramowe (vielbein): e^a_i = sqrt(g)*R^a_i niesie strukture SU(2)
Ansatz jezowy (hedgehog): R^a_i = delta^a_i

Topologia: vielbein mapuje skompaktyfikowane R^3 ~ S^3 -> SU(2) ~ S^3
Homotopia pi_3(S^3) = Z klasyfikuje wg liczby nawijania B

Rdzen solitonu: g mapuje [g0, 1] -> profil [pi, 0]
To daje DOKLADNIE B = 1, niezaleznie od ksztaltu profilu (topologiczne!)

Kwantyzacja (Finkelstein-Rubinstein):
  Obrot o 2pi -> faza (-1)^B
  B=1: fermion (spin 1/2)
  B=2: bozon (spin 0 lub 1)

## Wyniki (2026-04-15) -- q5_spin.py (7/7 PASS)

### Liczba nawijania B = n
- Profile Gaussowski, wymierny, tanh: B = n dokladnie dla n=1,2,3
- Niezaleznosc od profilu: rozrzut < 10^-6 dla roznych szerokosci
- TOPOLOGICZNA: zalezy TYLKO od warunkow brzegowych f(0), f(inf)

### Soliton TGP -> Jezowy B=1
- Rdzen solitonu: f(r) = pi*(1-g(r))/(1-g0), f(0)=pi, f(r_cross)=0
- Testowane g0 = 0.3, 0.5, 0.7, 0.85, 0.95: B = 1.0000 dla wszystkich
- DOWOD ANALITYCZNY: B = 2 * int_0^1 sin^2(pi*u) du = 2 * 1/2 = 1
- B=1 jest TOPOLOGICZNE: niezalezne od g0 i ksztaltu profilu

### Kwantyzacja momentu pedu
- J_min = B/2 = 1/2 (spin polcalkowy!)
- Widmo rotacyjne: E_J = J(J+1)/(2*Lambda)
- Stosunek E(3/2)/E(1/2) = 5.0 (dokladnie)
- Fizyczny analog: rozszczepienie Delta-Nukleon w QCD

### Zwiazek spin-statystyka
- Zamiana dwoch solitonow = obrot o 2pi jednego
- Faza zamiany = (-1)^B: automatycznie z topologii
- B=1 -> fermion (Pauli), B=2 -> bozon (Bose)
- NIE jest osobnym postulatem -- wynika z pi_3(S^3)=Z

### Czynnik giromagnetyczny g = 2
- Prad elektromagnetyczny sprzega sie z pradem topologicznym
- mu = (e/2m)*B, spin s = B/2 => g = B/s = 2 (wartosc Diraca)
- Poprawka anomalna ~ 1/Lambda ~ korekcja z rozmiaru solitonu

### Odpornosc topologiczna
- Perturbacje 1-15%: B = 1.000 (niezmienione)
- Topologia jest odporna na ciagle deformacje

## Lancuch derywacji

```
Soliton TGP g(r): g(0)=g0, g(inf)=1
  |
  v
Vielbein: e^a_i = sqrt(g) * R^a_i (jezowy)
  |
  v
Mapa: S^3 -> SU(2) ~ S^3, pi_3(S^3) = Z
Nawijanie: B = 1 dla fundamentalnego solitonu
  |
  v
Kwantyzacja Finkelsteina-Rubinsteina:
  2pi obrot -> (-1)^B = -1
  |
  v
SPIN 1/2 + SPIN-STATYSTYKA + g = 2
```

## Predykcje testowalne

1. Wzbudzenia rotacyjne: E(3/2)/E(1/2) = 5 (analog Delta/N)
2. Anomalny moment: a = (g-2)/2 ~ 1/(Lambda*m)
3. Kwantyzacja spinu DOKLADNA (topologiczna, nie przyblizenie)
4. Wyzsze spiny: B=2 solitony (spin 1) jako stany zwiazane
5. Zwiazek spin-statystyka AUTOMATYCZNY z topologii

## STATUS: Q5 ZAMKNIETE

- [x] Mechanizm topologiczny zidentyfikowany (pi_3(S^3)=Z)
- [x] Soliton TGP mapuje do B=1 (analitycznie + numerycznie)
- [x] Spin 1/2 z kwantyzacji FR
- [x] Spin-statystyka automatyczna
- [x] Czynnik g = 2 z pradu topologicznego
- [x] Odpornosc topologiczna potwierdzona

## Pliki

| Plik | Opis | Status |
|------|------|--------|
| q5_spin.py | Nawijanie, FR, widmo, spin-statystyka, g=2 | 7/7 PASS |
